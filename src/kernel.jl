# Banded pairwise alignment of encoded nucleotide sequences with affine gap
# penalties (Gotoh) and separate penalties for the gaps at the ends of each
# sequence (end gaps).
#
# The dynamic programming matrix is filled by anti-diagonals (`d = i + j`):
# the cells of an anti-diagonal depend only on the two previous ones, so the
# inner loop has no loop-carried dependency and is vectorised (SIMD). Scores
# are 16-bit integers when they cannot overflow (16 cells per AVX2
# instruction), 32-bit otherwise.
#
# Band coordinates: diagonal `k = j - i`, with `klo <= k <= khi`. On
# anti-diagonal `d` only diagonals with the parity of `d` occur; slot `m`
# (0-based) holds `k = kb(d) + 2m`, where `kb(d)` is the first such diagonal.

# Sequence `a` indexes the rows (`i`), `b` the columns (`j`). Horizontal gaps
# (gap symbols in `a`, consuming `b`) in the first or last row are end gaps of
# `a`; vertical gaps in the first or last column are end gaps of `b`.
#
# With `clip_open >= 0` (soft clipping) the first `i` bases of `a` may be left
# unaligned (the path starts at any cell of row `i`), and so may the last
# ones, at a cost of `clip_open + n * clip_extend` per clipped end. The end
# clip is a third vertical gap state G that can only return to H in the last
# row, from which the end gaps of `a` (free when `a` is fitted into `b`) lead
# to `(la, lb)`. Bases of `b` outside the path cost the end gaps of `a`.
struct _Scoring
    match::Int32
    mismatch::Int32
    gap_open::Int32
    gap_extend::Int32
    a_end_open::Int32
    a_end_extend::Int32
    b_end_open::Int32
    b_end_extend::Int32
    clip_open::Int32
    clip_extend::Int32
end

_Scoring(match, mismatch, gap_open, gap_extend, a_end_open, a_end_extend, b_end_open,
         b_end_extend) =
    _Scoring(match, mismatch, gap_open, gap_extend, a_end_open, a_end_extend, b_end_open,
             b_end_extend, -1, -1)

_clipping(sc::_Scoring) = sc.clip_open >= 0

struct _AlignStats
    score::Int
    matches::Int
    mismatches::Int
    columns::Int
    terminal_gap_columns::Int   # columns before the first / after the last aligned pair
    gap_openings::Int           # including end gaps (used by iddef 3)
    internal_gap_openings::Int  # between the first and the last aligned pair
    touches_band_edge::Bool     # the path reaches a band diagonal inside the matrix
    astart::Int                 # first and last positions of `a` and `b` in an
    aend::Int                   # aligned pair (0 if there is none)
    bstart::Int
    bend::Int
end

# Score buffers for one anti-diagonal: current (H, E, F, G), previous (H1,
# E1, F1, G1) and the one before (H2). G is used only with soft clipping.
mutable struct _Buffers{T<:Signed}
    H::Vector{T}; E::Vector{T}; F::Vector{T}; G::Vector{T}
    H1::Vector{T}; E1::Vector{T}; F1::Vector{T}; G1::Vector{T}
    H2::Vector{T}
end

_Buffers{T}() where {T} = _Buffers{T}(T[], T[], T[], T[], T[], T[], T[], T[], T[])

function _reset!(buf::_Buffers{T}, len::Int, neg::T) where {T}
    for v in (buf.H, buf.E, buf.F, buf.G, buf.H1, buf.E1, buf.F1, buf.G1, buf.H2)
        resize!(v, len)
        fill!(v, neg)
    end
    buf
end

# Per-task reusable memory.
mutable struct _Workspace
    buf16::_Buffers{Int16}
    buf32::_Buffers{Int32}
    trace::Vector{UInt8}
    ops::Vector{UInt8}
    arev::Vector{UInt8}   # `a` reversed, with a sentinel at each end
    bpad::Vector{UInt8}   # `b` with a sentinel at each end
end

_Workspace() = _Workspace(_Buffers{Int16}(), _Buffers{Int32}(), UInt8[], UInt8[], UInt8[], UInt8[])

# Int16 if no score (including the "minus infinity" sentinel minus a few
# penalties) can leave the Int16 range.
function _score_type(la::Int, lb::Int, sc::_Scoring)
    maxpen = max(abs(sc.mismatch), abs(sc.match), sc.gap_open + sc.gap_extend,
                 sc.a_end_open + sc.a_end_extend, sc.b_end_open + sc.b_end_extend)
    worst = abs(sc.mismatch) * min(la, lb) + sc.match * min(la, lb) +
            (sc.gap_open + max(sc.a_end_open, sc.b_end_open)) +
            max(sc.gap_extend, sc.a_end_extend, sc.b_end_extend) * (la + lb)
    _clipping(sc) && (worst += 2 * (sc.clip_open + sc.clip_extend * la);
                      maxpen = max(maxpen, sc.clip_open + sc.clip_extend))
    worst < 12_000 && maxpen < 1_000 ? Int16 : Int32
end

_kbase(d::Int, klo::Int) = klo + ((d - klo) & 1)

"""
    _align(ws, a, b, klo, khi, sc) -> _AlignStats

Alignment of encoded sequences `a` (rows `i`) and `b` (columns `j`) from
`(0, 0)` to `(la, lb)`, restricted to the diagonals `klo <= j - i <= khi`
(which must contain `0` and `length(b) - length(a)`). End gaps use the end
gap penalties of `sc`: with zero end gap penalties for `a`, `a` is fitted
into `b`. A gap of length `n` costs `open + n * extend`. Ties prefer the
diagonal, then a gap in `a`, then a gap in `b`; gap opening is preferred to
extension. With soft clipping (see `_Scoring`) a clipped end is used only if
it scores strictly better. The operations of the path are left in `ws.ops`,
in reverse order: `'M'` aligned pair, `'D'` gap in `a`, `'I'` gap in `b`
(clipped bases of `a` are `'I'`).
"""
_align(ws::_Workspace, a::AbstractVector{UInt8}, b::AbstractVector{UInt8}, klo::Int, khi::Int,
       sc::_Scoring) =
    _align(ws, a, b, klo, khi, sc, _score_type(length(a), length(b), sc))

function _align(ws::_Workspace, a::AbstractVector{UInt8}, b::AbstractVector{UInt8}, klo::Int,
                khi::Int, sc::_Scoring, ::Type{T}) where {T<:Signed}
    buf = T === Int16 ? ws.buf16 : ws.buf32
    score = _clipping(sc) ? _fill!(ws, buf::_Buffers{T}, a, b, klo, khi, sc, Val(true)) :
                            _fill!(ws, buf::_Buffers{T}, a, b, klo, khi, sc, Val(false))
    _traceback(ws, a, b, klo, khi, score)
end

# Fills the trace matrix: one byte per cell, bits 0-1 origin of H (0
# diagonal, 1 E, 2 F, 3 start after a clipped beginning of `a`), bit 2 E
# extends E, bit 3 F extends F, bit 4 G extends G, bit 5 H comes from G (last
# row only). Returns the score of the alignment.
function _fill!(ws::_Workspace, buf::_Buffers{T}, a::AbstractVector{UInt8},
                b::AbstractVector{UInt8}, klo::Int, khi::Int, sc::_Scoring,
                ::Val{C}) where {T,C}
    la = length(a); lb = length(b)
    M = (khi - klo) ÷ 2 + 1                  # slots per anti-diagonal
    neg = typemin(T) ÷ T(2)
    _reset!(buf, M + 2, neg)                 # slot m at index m + 2; padding at both ends
    ncell = (la + lb + 1) * M
    length(ws.trace) < ncell && resize!(ws.trace, ncell)
    trace = ws.trace
    arev = resize!(ws.arev, la + 2)          # a[i] at index la - i + 2
    arev[1] = 0x05; arev[la+2] = 0x05
    for i in 1:la
        arev[la-i+2] = a[i]
    end
    bpad = resize!(ws.bpad, lb + 2)          # b[j] at index j + 1
    bpad[1] = 0x05; bpad[lb+2] = 0x05
    copyto!(bpad, 2, b, 1, lb)

    match = T(sc.match); mismatch = T(sc.mismatch)
    gap_oe = T(sc.gap_open + sc.gap_extend); gap_e = T(sc.gap_extend)
    co = T(max(sc.clip_open, 0)); ce = T(max(sc.clip_extend, 0)); clip_oe = co + ce
    H, E, F, G, H1, E1, F1, G1, H2 = buf.H, buf.E, buf.F, buf.G, buf.H1, buf.E1, buf.F1,
                                     buf.G1, buf.H2
    @inbounds for d in 0:(la + lb)
        # Rotate: d-1 becomes "previous", d-2 "before previous".
        H, H1, H2 = H2, H, H1
        E, E1 = E1, E
        F, F1 = F1, F
        G, G1 = G1, G
        kb = _kbase(d, klo)
        s = kb == klo ? -1 : 0               # slot shift of the neighbours on d - 1
        kmin = max(klo, -d, d - 2la)
        kmax = min(khi, d, 2lb - d)
        mlo = (kmin - kb + 1) >> 1           # cld(kmin - kb, 2)
        mhi = (kmax - kb) >> 1               # fld(kmax - kb, 2)
        # Slots outside the matrix hold minus infinity.
        for m in 0:min(mlo - 1, M - 1)
            H[m+2] = neg; E[m+2] = neg; F[m+2] = neg
            C && (G[m+2] = neg)
        end
        for m in max(mhi + 1, 0):(M - 1)
            H[m+2] = neg; E[m+2] = neg; F[m+2] = neg
            C && (G[m+2] = neg)
        end
        mlo > mhi && continue
        ib = (d - kb) >> 1                   # row i of slot m is ib - m
        abase = la - ib + 2
        bbase = ((d + kb) >> 1) + 1
        tbase = d * M + 1
        # Interior cells with internal gap penalties (vectorised).
        @simd ivdep for m in mlo:mhi
            x = m + 2
            e_open = H1[x+s] - gap_oe
            e_ext = E1[x+s] - gap_e
            eext = e_ext > e_open
            e = ifelse(eext, e_ext, e_open)
            f_open = H1[x+s+1] - gap_oe
            f_ext = F1[x+s+1] - gap_e
            fext = f_ext > f_open
            f = ifelse(fext, f_ext, f_open)
            ai = arev[abase+m]
            bj = bpad[bbase+m]
            diag = H2[x] + ifelse((ai == bj) & (ai < 0x04), match, mismatch)
            src = ifelse(e > diag, 0x01, 0x00)
            h = max(diag, e)
            src = ifelse(f > h, 0x02, src)
            hf = max(h, f)
            tb = (UInt8(eext) << 2) | (UInt8(fext) << 3)
            if C
                # Start after clipping the first i bases of `a` (`% T`: no
                # checked conversion, which would prevent vectorisation).
                cs = -(co + ce * ((ib - m) % T))
                src = ifelse(cs > hf, 0x03, src)
                hf = max(hf, cs)
                # End clip state (vertical, from the cell above).
                g_open = H1[x+s+1] - clip_oe
                g_ext = G1[x+s+1] - ce
                gext = g_ext > g_open
                G[x] = ifelse(gext, g_ext, g_open)
                tb |= UInt8(gext) << 4
            end
            H[x] = hf; E[x] = e; F[x] = f
            trace[tbase+m] = src | tb
        end
        # Cells in the first/last row or column (end gaps, origin) lie at the
        # ends of the slot range: recompute them.
        _boundary_cell!(H, E, F, G, H1, E1, F1, G1, H2, trace, arev, bpad, d, kb, mlo, s, la,
                        lb, abase, bbase, tbase, neg, sc, Val(C), co, ce)
        mhi > mlo && _boundary_cell!(H, E, F, G, H1, E1, F1, G1, H2, trace, arev, bpad, d, kb,
                                      mhi, s, la, lb, abase, bbase, tbase, neg, sc, Val(C), co,
                                      ce)
    end
    # After the last rotation H holds anti-diagonal la + lb.
    Int(H[((lb - la - _kbase(la + lb, klo)) >> 1) + 2])
end

@inline function _boundary_cell!(H::Vector{T}, E, F, G, H1, E1, F1, G1, H2, trace, arev, bpad,
                                 d, kb, m, s, la, lb, abase, bbase, tbase, neg::T,
                                 sc::_Scoring, ::Val{C}, co::T, ce::T) where {T,C}
    k = kb + 2m
    i = (d - k) >> 1
    j = (d + k) >> 1
    rowterm = i == 0 || i == la
    colterm = j == 0 || j == lb
    (rowterm || colterm) || return nothing
    x = m + 2
    @inbounds if i == 0 && j == 0
        H[x] = 0; E[x] = neg; F[x] = neg
        C && (G[x] = neg)
        trace[tbase+m] = 0x00
        return nothing
    end
    hgo = T(rowterm ? sc.a_end_open : sc.gap_open)
    hge = T(rowterm ? sc.a_end_extend : sc.gap_extend)
    vgo = T(colterm ? sc.b_end_open : sc.gap_open)
    vge = T(colterm ? sc.b_end_extend : sc.gap_extend)
    @inbounds begin
        e_open = H1[x+s] - hgo - hge
        e_ext = E1[x+s] - hge
        eext = e_ext > e_open
        e = ifelse(eext, e_ext, e_open)
        f_open = H1[x+s+1] - vgo - vge
        f_ext = F1[x+s+1] - vge
        fext = f_ext > f_open
        f = ifelse(fext, f_ext, f_open)
        ai = arev[abase+m]
        bj = bpad[bbase+m]
        diag = H2[x] + T(((ai == bj) & (ai < 0x04)) ? sc.match : sc.mismatch)
        src = ifelse(e > diag, 0x01, 0x00)
        h = max(diag, e)
        src = ifelse(f > h, 0x02, src)
        hf = max(h, f)
        tb = (UInt8(eext) << 2) | (UInt8(fext) << 3)
        if C
            if i >= 1
                cs = -(co + ce * T(i))
                src = ifelse(cs > hf, 0x03, src)
                hf = max(hf, cs)
                g_open = H1[x+s+1] - co - ce
                g_ext = G1[x+s+1] - ce
                gext = g_ext > g_open
                g = ifelse(gext, g_ext, g_open)
                tb |= UInt8(gext) << 4
                if i == la && g > hf         # end clip returns to H in the last row
                    hf = g
                    tb |= 0x20
                end
                G[x] = g
            else
                G[x] = neg
            end
        end
        H[x] = hf; E[x] = e; F[x] = f
        trace[tbase+m] = src | tb
    end
    nothing
end

# Follows the trace from (la, lb) to (0, 0); operations: 'M' aligned pair,
# 'D' gap in `a` (consumes b), 'I' gap in `b` (consumes a), stored in reverse.
# Clipped bases of `a` are 'I'; the bases of `b` before a clipped start are
# 'D'.
function _traceback(ws::_Workspace, a::AbstractVector{UInt8}, b::AbstractVector{UInt8}, klo::Int,
                    khi::Int, score::Int)
    la = length(a); lb = length(b)
    M = (khi - klo) ÷ 2 + 1
    trace = ws.trace
    ops = empty!(ws.ops)
    i = la; j = lb; state = 0
    matches = 0
    edge = false
    astart = aend = bstart = bend = 0
    @inbounds while i > 0 || j > 0
        k = j - i
        d = i + j
        edge |= (k == klo && klo > -la) | (k == khi && khi < lb)
        tb = trace[d*M + ((k - _kbase(d, klo)) >> 1) + 1]
        if state == 0
            src = tb & 0x03
            if (tb & 0x20) != 0              # H from the end clip state
                state = 3
            elseif src == 0x03               # start after clipping the first i bases of a
                for _ in 1:i; push!(ops, UInt8('I')); end
                for _ in 1:j; push!(ops, UInt8('D')); end
                i = 0; j = 0
            elseif src == 0x00
                (a[i] == b[j] && a[i] < 0x04) && (matches += 1)
                aend == 0 && (aend = i; bend = j)
                astart = i; bstart = j
                push!(ops, UInt8('M')); i -= 1; j -= 1
            else
                state = src == 0x01 ? 1 : 2
            end
        elseif state == 1
            push!(ops, UInt8('D'))
            state = (tb & 0x04) != 0 ? 1 : 0
            j -= 1
        elseif state == 2
            push!(ops, UInt8('I'))
            state = (tb & 0x08) != 0 ? 2 : 0
            i -= 1
        else
            push!(ops, UInt8('I'))
            state = (tb & 0x10) != 0 ? 3 : 0
            i -= 1
        end
    end

    ncol = length(ops)
    first_m = findfirst(==(UInt8('M')), ops)
    aligned = count(==(UInt8('M')), ops)
    terminal = first_m === nothing ? ncol :
               (first_m - 1) + (ncol - findlast(==(UInt8('M')), ops))
    openings = 0; internal = 0
    lo = first_m === nothing ? 1 : first_m
    hi = first_m === nothing ? 0 : findlast(==(UInt8('M')), ops)
    @inbounds for k in 1:ncol
        o = ops[k]
        if o != UInt8('M') && (k == 1 || ops[k-1] != o)
            openings += 1
            lo < k < hi && (internal += 1)
        end
    end
    _AlignStats(score, matches, aligned - matches, ncol, terminal, openings, internal, edge,
                astart, aend, bstart, bend)
end

# CIGAR string of the path in `ops` (reverse order) from the first to the last
# aligned pair: `M` aligned pair, `I` base of `a` only, `D` base of `b` only;
# unaligned bases of `a` before / after are soft clips (`S`); unaligned bases
# of `b` are left out (they are given by the coordinates).
function _cigar(ops::Vector{UInt8}, st::_AlignStats, la::Int)
    st.aend == 0 && return string(la, 'S')
    io = IOBuffer()
    st.astart > 1 && print(io, st.astart - 1, 'S')
    lo = findfirst(==(UInt8('M')), ops); hi = findlast(==(UInt8('M')), ops)
    k = hi
    while k >= lo
        o = ops[k]; n = 0
        while k >= lo && ops[k] == o
            n += 1; k -= 1
        end
        print(io, n, Char(o))
    end
    st.aend < la && print(io, la - st.aend, 'S')
    String(take!(io))
end
