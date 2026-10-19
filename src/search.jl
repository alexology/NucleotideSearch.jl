# Search strategy of `vsearch --usearch_global`: candidate targets are ranked
# by the number of unique k-mers shared with the query and aligned in that
# order until `maxaccepts` hits are accepted or `maxrejects` candidates are
# rejected.

"""
    Hit

One accepted alignment of a search. Fields: `query`, `target` (indices),
`strand` (`'+'`, or `'-'` if the reverse complement of the query was
aligned), `identity` (percent, see `iddef` in [`search`](@ref)), `score`,
`matches`, `mismatches`, `gap_openings`, `columns`, `aligned_columns`,
`query_start`, `query_end`, `target_start`, `target_end`, `cigar`.

`aligned_columns` counts the alignment columns and `gap_openings` the gaps
from the first to the last aligned pair (end gaps excluded; `iddef = 3`
still counts end gaps, as VSEARCH). `columns` also counts the end gaps; in `:query` mode the
unaligned target bases counted there depend on the window used, so use
`aligned_columns` (plus unaligned query bases, if wanted) instead.

Coordinates are 1-based and inclusive and give the first and last bases in an
aligned pair. Query coordinates refer to the query as given (also on the
`'-'` strand); target coordinates always satisfy `target_start <=
target_end`. The CIGAR string (empty unless requested) describes the aligned
strand of the query against the target, from `target_start` to
`target_end`: `M` aligned pair, `I` query base only, `D` target base only,
`S` unaligned query bases at the ends.
"""
struct Hit
    query::Int
    target::Int
    strand::Char
    identity::Float64
    score::Int
    matches::Int
    mismatches::Int
    gap_openings::Int
    columns::Int
    aligned_columns::Int
    query_start::Int
    query_end::Int
    target_start::Int
    target_end::Int
    cigar::String
end

function _identity(st::_AlignStats, la::Int, lb::Int, iddef::Int)
    if iddef == 0
        return 100 * st.matches / min(la, lb)
    elseif iddef == 1 || iddef == 4
        return 100 * st.matches / st.columns
    elseif iddef == 2
        den = st.columns - st.terminal_gap_columns
        return den == 0 ? 0.0 : 100 * st.matches / den
    else
        return 100 * (1 - (st.mismatches + st.gap_openings) / max(la, lb))
    end
end

# Settings shared by all queries of a search.
struct _SearchSettings
    mode::Symbol
    min_identity::Float64   # fraction
    id_threshold::Float64   # percent, with a tolerance for rounding
    iddef::Int
    query_cov::Float64
    maxhits::Int
    maxaccepts::Int
    maxrejects::Int
    minwordmatches::Int
    scoring::_Scoring
    cigar::Bool
end

function _settings(mode, min_identity, iddef, query_cov, maxhits, maxaccepts, maxrejects,
                   minwordmatches, scoring::Scoring, cigar)
    0 <= min_identity <= 100 ||
        throw(ArgumentError("min_identity must be a percentage between 0 and 100"))
    0 < min_identity <= 1 &&
        @warn "min_identity = $min_identity is a percentage: did you mean $(100 * min_identity)?"
    iddef in 0:4 || throw(ArgumentError("iddef must be 0, 1, 2, 3 or 4"))
    0 <= query_cov <= 1 || throw(ArgumentError("query_cov must be between 0 and 1"))
    min(maxhits, maxaccepts, maxrejects, minwordmatches) >= 0 ||
        throw(ArgumentError("maxhits, maxaccepts, maxrejects and minwordmatches must be >= 0"))
    _SearchSettings(mode, min_identity / 100, Float64(min_identity) - 1e-9, iddef, query_cov,
                    maxhits, maxaccepts, maxrejects, minwordmatches,
                    _kernel_scoring(scoring, mode), cigar)
end

# Diagonal band for sequences of lengths `la` and `lb`. Under iddef 1/4 an
# alignment with identity >= t has at most `slack` gap columns, so it stays
# within (slack -/+ (lb - la)) / 2 of the main diagonal. Other definitions get a
# wider band.
function _band(la::Int, lb::Int, s::_SearchSettings)
    t = s.min_identity
    slack = t > 0 ? floor(Int, (1 - t) / t * min(la, lb)) : la + lb
    delta = lb - la
    if s.iddef == 1 || s.iddef == 4
        max(-la, min(0, delta, -fld(slack - delta, 2))), min(lb, max(0, delta, fld(slack + delta, 2)))
    else
        max(-la, min(0, delta) - slack), min(lb, max(0, delta) + slack)
    end
end

# Exact test on the lengths alone: false if no alignment can be accepted.
function _length_compatible(la::Int, lb::Int, s::_SearchSettings)
    la > 0 && lb > 0 && min(la, lb) >= s.query_cov * la || return false
    s.mode === :query && return true
    !((s.iddef == 1 || s.iddef == 4) && 100 * min(la, lb) / max(la, lb) < s.id_threshold)
end

# Per-task memory of a search.
struct _SearchState
    ws::_Workspace
    shared::Vector{Int32}    # shared unique k-mers with the current query
    touched::Vector{Int}
    candidates::Vector{Int}
    qkp::Vector{UInt32}      # k-mers occurring once in the query, sorted,
    qpos::Vector{Int32}      # with their end positions
    diags::Vector{Int}
end

_SearchState(n::Int) = _SearchState(_Workspace(), zeros(Int32, n), Int[], Int[], UInt32[],
                                    Int32[], Int[])

# Targets sharing enough unique k-mers (`qkmers`) with the query, by decreasing
# number of shared k-mers (ties: index).
function _rank_candidates!(st::_SearchState, qkmers::Vector{UInt32}, db::SequenceDB,
                           s::_SearchSettings)
    shared = st.shared; touched = empty!(st.touched); candidates = empty!(st.candidates)
    offsets = db.offsets; postings = db.postings
    for w in qkmers, p in (offsets[w+1]+1):offsets[w+2]
        tg = postings[p]
        shared[tg] == 0 && push!(touched, tg)
        shared[tg] += Int32(1)
    end
    threshold = min(s.minwordmatches, length(qkmers))
    if threshold == 0
        append!(candidates, eachindex(shared))
    else
        append!(candidates, Iterators.filter(tg -> shared[tg] >= threshold, touched))
    end
    sort!(candidates; alg = MergeSort, by = tg -> (-shared[tg], tg))
    for tg in touched
        shared[tg] = 0
    end
    candidates
end

function _hit(q::Int, tg::Int, strand::Char, id::Float64, st::_AlignStats, la::Int, cigar::String)
    qs, qe = strand == '+' ? (st.astart, st.aend) :
             st.aend == 0 ? (0, 0) : (la - st.aend + 1, la - st.astart + 1)
    Hit(q, tg, strand, id, st.score, st.matches, st.mismatches, st.internal_gap_openings, st.columns,
        st.columns - st.terminal_gap_columns, qs, qe, st.bstart, st.bend, cigar)
end

# Hits are sorted by decreasing identity (stable: candidate order on ties).
function _finish!(hits::Vector{Hit}, s::_SearchSettings)
    sort!(hits; alg = MergeSort, by = h -> -h.identity)
    s.maxhits > 0 && length(hits) > s.maxhits && resize!(hits, s.maxhits)
    hits
end

function _run_parallel(f, n::Int, ntargets::Int)
    chunks = Iterators.partition(1:n, max(1, cld(n, 8 * Threads.nthreads())))
    tasks = map(chunks) do range
        Threads.@spawn begin
            state = _SearchState(ntargets)
            for q in range
                f(state, q)
            end
        end
    end
    foreach(wait, tasks)
end

# ---------------------------------------------------------------------------
# All against all

# Alignments shared between the two directions of a pair, in lock-striped
# dictionaries. An entry is removed when it is used the second time.
struct _PairCache
    tables::Vector{Dict{Int,_AlignStats}}
    locks::Vector{ReentrantLock}
end

_PairCache(nstripes::Int) =
    _PairCache([Dict{Int,_AlignStats}() for _ in 1:nstripes], [ReentrantLock() for _ in 1:nstripes])

function _get_or_align!(compute, cache::_PairCache, key::Int)
    stripe = mod1(key, length(cache.tables))
    table = cache.tables[stripe]
    lk = cache.locks[stripe]
    cached = lock(() -> pop!(table, key, nothing), lk)
    cached === nothing || return cached
    st = compute()   # outside the lock; a concurrent duplicate gives the same result
    lock(lk) do
        haskey(table, key) ? delete!(table, key) : (table[key] = st)
    end
    st
end

# Alignment of sequences `x` < `y`. A pair is always aligned in this canonical
# orientation, so (q, t) and (t, q) share one alignment (identity and number of
# aligned pairs are symmetric). An accepted alignment that touches the band edge
# may not be optimal: it is recomputed without band.
function _pair_stats(ws::_Workspace, db::SequenceDB, x::Int, y::Int, s::_SearchSettings)
    a = _seq(db, x); b = _seq(db, y)
    la = length(a); lb = length(b)
    klo, khi = _band(la, lb, s)
    st = _align(ws, a, b, klo, khi, s.scoring)
    if st.touches_band_edge && _identity(st, la, lb, s.iddef) >= s.id_threshold
        st = _align(ws, a, b, -la, lb, s.scoring)
    end
    st
end

# The pair was aligned as (x, y) = (t, q): swap the roles back.
_swap(st::_AlignStats) =
    _AlignStats(st.score, st.matches, st.mismatches, st.columns, st.terminal_gap_columns,
                st.gap_openings, st.internal_gap_openings, st.touches_band_edge, st.bstart,
                st.bend, st.astart, st.aend)

function _search_self_query(st::_SearchState, q::Int, db::SequenceDB, cache::_PairCache,
                            s::_SearchSettings)
    n = nsequences(db)
    la = _seqlen(db, q)
    accepted = Hit[]
    rejected = 0
    for tg in _rank_candidates!(st, _unique_kmers(_seq(db, q), db.k), db, s)
        tg == q && continue
        lb = _seqlen(db, tg)
        ok = false
        if _length_compatible(la, lb, s)
            x, y = minmax(q, tg)
            stats = _get_or_align!(cache, (x - 1) * n + y) do
                _pair_stats(st.ws, db, x, y, s)
            end
            id = _identity(stats, la, lb, s.iddef)
            ok = id >= s.id_threshold && stats.matches + stats.mismatches >= s.query_cov * la
            ok && push!(accepted, _hit(q, tg, '+', id, q < tg ? stats : _swap(stats), la, ""))
        end
        if ok
            s.maxaccepts > 0 && length(accepted) >= s.maxaccepts && break
        else
            rejected += 1
            s.maxrejects > 0 && rejected >= s.maxrejects && break
        end
    end
    _finish!(accepted, s)
end

"""
    search_self(db::SequenceDB; min_identity = 84, iddef = 1, query_cov = 0.9,
                maxhits = 10, maxaccepts = 0, maxrejects = 32,
                minwordmatches = 12, scoring = Scoring()) -> Vector{Vector{Hit}}

Search every sequence of `db` (query) against all the others (targets) with
a global alignment, as `vsearch --usearch_global --self`; element `q` of the
result holds the hits of sequence `q`. Self-hits are never reported.

For each query:

1. targets are ranked by the number of unique k-mers (`db.k`) shared with the
   query, decreasing; targets sharing fewer than `minwordmatches` (or all of
   the query's unique k-mers, if fewer) are ignored;
2. candidates are aligned in that order with a banded global alignment (see
   [`Scoring`](@ref));
3. a candidate is accepted if its identity (percent) is `>= min_identity` and
   the fraction of query bases aligned to target bases is `>= query_cov`;
   otherwise it is rejected. The search stops after `maxaccepts` accepts (0:
   no limit) or `maxrejects` rejects (0: no limit);
4. the accepted hits are sorted by decreasing identity and the first
   `maxhits` (0: all) are kept.

With `maxrejects = 0` and `minwordmatches = 0` every pair is aligned
(exhaustive search; slow for large data sets).

Identity definitions (`iddef`, as in VSEARCH): `0`: matches / length of the
shorter sequence; `1` and `4`: matches / alignment columns; `2`: matches /
alignment columns excluding end gaps; `3`: `1 - (mismatches + gap openings) /
length of the longer sequence`. Bases other than A, C, G, T/U never count as
matches.

The alignment is restricted to a diagonal band wide enough for every
alignment with identity `>= min_identity` under `iddef` 1. Queries run in
parallel on `Threads.nthreads()` threads; the result does not depend on the
number of threads.
"""
function search_self(db::SequenceDB; min_identity::Real = 84, iddef::Integer = 1,
                     query_cov::Real = 0.9, maxhits::Integer = 10, maxaccepts::Integer = 0,
                     maxrejects::Integer = 32, minwordmatches::Integer = 12,
                     scoring::Scoring = Scoring())
    s = _settings(:global, min_identity, iddef, query_cov, maxhits, maxaccepts, maxrejects,
                  minwordmatches, scoring, false)
    n = nsequences(db)
    cache = _PairCache(64)
    results = Vector{Vector{Hit}}(undef, n)
    _run_parallel(n, n) do state, q
        results[q] = _search_self_query(state, q, db, cache, s)
    end
    results
end

# ---------------------------------------------------------------------------
# Queries against a database

# Most frequent diagonal (target position - query position) among the k-mers
# that occur once in the query and appear in the target; `nothing` if none.
function _seed_diagonal(st::_SearchState, t::AbstractVector{UInt8}, k::Int)
    kp = st.qkp; pos = st.qpos
    isempty(kp) && return nothing
    diags = empty!(st.diags)
    mask = (UInt32(1) << (2k)) - UInt32(1)
    w = UInt32(0); valid = 0
    @inbounds for (p, x) in enumerate(t)
        if x > 0x03
            valid = 0; w = UInt32(0)
            continue
        end
        w = ((w << 2) | x) & mask
        valid += 1
        if valid >= k
            r = searchsortedfirst(kp, w)
            r <= length(kp) && kp[r] == w && push!(diags, p - Int(pos[r]))
        end
    end
    isempty(diags) && return nothing
    sort!(diags)
    best = diags[1]; bestn = 0; i = 1
    @inbounds while i <= length(diags)
        j = i
        while j < length(diags) && diags[j+1] == diags[i]
            j += 1
        end
        j - i + 1 > bestn && (best = diags[i]; bestn = j - i + 1)
        i = j + 1
    end
    best
end

# Fitting (`:query` mode) of query `a` into target `b`, in a window of the
# target around the seed diagonal. The result is recomputed on the whole target
# without band if the path touches the band edge or a window boundary inside
# the target. Coordinates refer to the whole target.
function _query_stats(st::_SearchState, a::AbstractVector{UInt8}, b::AbstractVector{UInt8},
                      k::Int, s::_SearchSettings)
    la = length(a); lb = length(b)
    t = s.min_identity
    slack = t > 0 ? min(la + lb, floor(Int, (1 - t) / t * la)) : la + lb
    d = _seed_diagonal(st, b, k)
    if d !== nothing
        margin = slack + 16
        lo = max(1, d + 1 - margin); hi = min(lb, d + la + margin)
        if lo <= hi
            w = view(b, lo:hi); lw = length(w)
            klo = max(-la, -slack); khi = min(lw, lw - la + slack)
            if klo <= min(0, lw - la) && khi >= max(0, lw - la)
                r = _align(st.ws, a, w, klo, khi, s.scoring)
                inside = r.aend > 0 && (lo == 1 || r.bstart > 1) && (hi == lb || r.bend < lw)
                if inside && !r.touches_band_edge
                    return _AlignStats(r.score, r.matches, r.mismatches, r.columns,
                                       r.terminal_gap_columns, r.gap_openings,
                                       r.internal_gap_openings, false,
                                       r.astart, r.aend, r.bstart + lo - 1, r.bend + lo - 1)
                end
            end
        end
    end
    _align(st.ws, a, b, -la, lb, s.scoring)
end

function _target_stats(st::_SearchState, a, b, k::Int, s::_SearchSettings)
    if s.mode === :query
        _query_stats(st, a, b, k, s)
    else
        la = length(a); lb = length(b)
        klo, khi = _band(la, lb, s)
        r = _align(st.ws, a, b, klo, khi, s.scoring)
        if r.touches_band_edge && _identity(r, la, lb, s.iddef) >= s.id_threshold
            r = _align(st.ws, a, b, -la, lb, s.scoring)
        end
        r
    end
end

function _search_query(st::_SearchState, q::Int, query::Vector{UInt8}, db::SequenceDB,
                       strands::Tuple, s::_SearchSettings)
    la = length(query)
    k = db.k
    oriented = Tuple(strand == '+' ? query : reverse_complement(query) for strand in strands)
    # Candidates of all strands, by decreasing shared k-mers (ties: target, strand).
    cands = Tuple{Int32,Int,Int}[]   # (-shared, target, strand index)
    for (si, a) in enumerate(oriented)
        _collect_candidates!(cands, st, _unique_kmers(a, k), db, s, si)
    end
    sort!(cands)
    accepted = Hit[]
    rejected = 0
    current = 0
    for (_, tg, si) in cands
        b = _seq(db, tg)
        lb = length(b)
        a = oriented[si]
        ok = false
        if _length_compatible(la, lb, s)
            if s.mode === :query && si != current
                _kmer_positions!(st.qkp, st.qpos, a, k); current = si
            end
            r = _target_stats(st, a, b, k, s)
            id = _identity(r, la, lb, s.iddef)
            ok = r.aend > 0 && id >= s.id_threshold && r.matches + r.mismatches >= s.query_cov * la
            # The workspace still holds the path of `r`.
            ok && push!(accepted, _hit(q, tg, strands[si], id, r, la,
                                       s.cigar ? _cigar(st.ws.ops, r, la) : ""))
        end
        if ok
            s.maxaccepts > 0 && length(accepted) >= s.maxaccepts && break
        else
            rejected += 1
            s.maxrejects > 0 && rejected >= s.maxrejects && break
        end
    end
    _finish!(accepted, s)
end

# Appends (-shared k-mers, target, `si`) for the targets sharing enough of the
# query k-mers `qk` (all targets if the threshold is 0).
function _collect_candidates!(cands, st::_SearchState, qk::Vector{UInt32}, db::SequenceDB,
                              s::_SearchSettings, si::Int)
    shared = st.shared; touched = empty!(st.touched)
    offsets = db.offsets; postings = db.postings
    for w in qk, p in (offsets[w+1]+1):offsets[w+2]
        tg = postings[p]
        shared[tg] == 0 && push!(touched, tg)
        shared[tg] += Int32(1)
    end
    threshold = min(s.minwordmatches, length(qk))
    if threshold == 0
        for tg in eachindex(shared)
            push!(cands, (-shared[tg], tg, si))
        end
    else
        for tg in touched
            shared[tg] >= threshold && push!(cands, (-shared[tg], tg, si))
        end
    end
    for tg in touched
        shared[tg] = 0
    end
    nothing
end

"""
    search(queries, db::SequenceDB; mode = :query, strand = :both,
           min_identity = 80, iddef = 2, query_cov = 0.9, maxhits = 10,
           maxaccepts = 10, maxrejects = 32, minwordmatches = 12,
           scoring = <mode default>, cigar = false) -> Vector{Vector{Hit}}

Search each query (strings or encoded vectors) against the database; element
`q` of the result holds the [`Hit`](@ref)s of query `q`, sorted by
decreasing identity.

- `mode = :query` (default): the whole query is fitted into the target; gaps
  at the ends of the query (unaligned target bases) are free, so target
  coordinates are exact start and end positions of the query on the
  reference, also with mismatches or indels near the ends. In the default
  scoring unaligned query bases beyond an end of the target cost
  `terminal_gap_open = 20` + `terminal_gap_extend = 1` per base, and so does
  a query end left unaligned inside the target (soft clipping): an end is
  clipped after about 7 consecutive mismatches (about 14 unrelated bases),
  so a divergent end does not distort the coordinates of the rest. Only
  `iddef` 0 and 2 are allowed. The alignment is computed in a window
  of the target around the most frequent k-mer diagonal and recomputed on the
  whole target if it reaches the window or band boundaries.
- `mode = :global`: end-to-end alignment as in [`search_self`](@ref).
- `strand = :both` also searches the reverse complement of each query
  (`:plus`: forward only).

Candidate selection, acceptance (`min_identity`, `query_cov`), `maxaccepts`,
`maxrejects`, `maxhits` and `iddef` are as in [`search_self`](@ref); `iddef =
2` (default here) ignores the unaligned ends of the target. Candidates of both
strands are ranked together. `cigar = true` adds the CIGAR string of each hit.
"""
function search(queries::AbstractVector, db::SequenceDB; mode::Symbol = :query,
                strand::Symbol = :both, min_identity::Real = 80, iddef::Integer = 2,
                query_cov::Real = 0.9, maxhits::Integer = 10, maxaccepts::Integer = 10,
                maxrejects::Integer = 32, minwordmatches::Integer = 12,
                scoring::Scoring = _default_scoring(_check_mode(mode)), cigar::Bool = false)
    _check_mode(mode)
    strand in (:both, :plus) || throw(ArgumentError("strand must be :both or :plus"))
    # Columns and gap openings that include end gaps depend on the target window
    # in `:query` mode, so only iddef 0 and 2 are meaningful there.
    mode === :query && !(iddef in (0, 2)) &&
        throw(ArgumentError("in :query mode iddef must be 0 or 2"))
    strands = strand === :both ? ('+', '-') : ('+',)
    s = _settings(mode, min_identity, iddef, query_cov, maxhits, maxaccepts, maxrejects,
                  minwordmatches, scoring, cigar)
    encoded = [collect(encode(x)) for x in queries]
    n = length(encoded)
    results = Vector{Vector{Hit}}(undef, n)
    _run_parallel(n, nsequences(db)) do state, q
        results[q] = _search_query(state, q, encoded[q], db, strands, s)
    end
    results
end

search(query::Union{AbstractString,AbstractVector{UInt8}}, db::SequenceDB; kwargs...) =
    only(search([query], db; kwargs...))
