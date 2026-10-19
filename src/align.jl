"""
    Scoring(; match = 2, mismatch = -4, gap_open = 20, gap_extend = 2,
            terminal_gap_open = 2, terminal_gap_extend = 1,
            query_end_gap_open = terminal_gap_open,
            query_end_gap_extend = terminal_gap_extend,
            target_end_gap_open = terminal_gap_open,
            target_end_gap_extend = terminal_gap_extend)

Alignment scores. Gap penalties are non-negative numbers; a gap of length
`n` costs `open + n * extend`. End gaps are gaps before the first or after
the last base of a sequence: `query_end_gap_*` apply to gap symbols at the
ends of the query (unaligned target bases), `target_end_gap_*` to gap symbols
at the ends of the target (unaligned query bases). The defaults are those of
VSEARCH.

In `:query` mode the query end gaps are always free, whatever the values
given here.
"""
struct Scoring
    match::Int
    mismatch::Int
    gap_open::Int
    gap_extend::Int
    query_end_gap_open::Int
    query_end_gap_extend::Int
    target_end_gap_open::Int
    target_end_gap_extend::Int
    function Scoring(match, mismatch, gap_open, gap_extend, query_end_gap_open,
                     query_end_gap_extend, target_end_gap_open, target_end_gap_extend)
        min(gap_open, gap_extend, query_end_gap_open, query_end_gap_extend,
            target_end_gap_open, target_end_gap_extend) >= 0 ||
            throw(ArgumentError("gap penalties must be given as non-negative numbers"))
        max(abs(match), abs(mismatch), gap_open, gap_extend, query_end_gap_open,
            query_end_gap_extend, target_end_gap_open, target_end_gap_extend) <= 10_000 ||
            throw(ArgumentError("scores and penalties must not exceed 10,000 in absolute value"))
        new(match, mismatch, gap_open, gap_extend, query_end_gap_open, query_end_gap_extend,
            target_end_gap_open, target_end_gap_extend)
    end
end

Scoring(; match::Integer = 2, mismatch::Integer = -4, gap_open::Integer = 20,
        gap_extend::Integer = 2, terminal_gap_open::Integer = 2,
        terminal_gap_extend::Integer = 1,
        query_end_gap_open::Integer = terminal_gap_open,
        query_end_gap_extend::Integer = terminal_gap_extend,
        target_end_gap_open::Integer = terminal_gap_open,
        target_end_gap_extend::Integer = terminal_gap_extend) =
    Scoring(match, mismatch, gap_open, gap_extend, query_end_gap_open, query_end_gap_extend,
            target_end_gap_open, target_end_gap_extend)

# Default scoring of each mode. In `:query` mode unaligned query bases at the
# ends (target too short, or ends that do not match it) cost almost as much as
# an internal gap (open 20, extend 1 instead of 2): the ends of the query are
# clipped only after about 8 consecutive mismatching bases, and an overhang
# beyond the end of the target is an end gap rather than an internal gap.
_default_scoring(mode::Symbol) =
    mode === :query ? Scoring(terminal_gap_open = 20, terminal_gap_extend = 1) : Scoring()

function _check_mode(mode::Symbol)
    mode in (:global, :query) || throw(ArgumentError("mode must be :global or :query"))
    mode
end

# Kernel scoring with the query as `a` (rows) and the target as `b`.
_kernel_scoring(sc::Scoring, mode::Symbol) =
    _Scoring(sc.match, sc.mismatch, sc.gap_open, sc.gap_extend,
             mode === :query ? 0 : sc.query_end_gap_open,
             mode === :query ? 0 : sc.query_end_gap_extend,
             sc.target_end_gap_open, sc.target_end_gap_extend)

"""
    PairAlignment

Result of [`align`](@ref): `score`, `matches`, `mismatches`, `gap_openings`
(between the first and the last aligned pair), `columns` (alignment
columns, including end gaps), `query_start`,
`query_end`, `target_start`, `target_end` (1-based, inclusive, first and last
bases in an aligned pair; `0` if there is none), `cigar` (see [`Hit`](@ref)).
"""
struct PairAlignment
    score::Int
    matches::Int
    mismatches::Int
    gap_openings::Int
    columns::Int
    query_start::Int
    query_end::Int
    target_start::Int
    target_end::Int
    cigar::String
end

"""
    align(query, target; mode = :global, scoring = <mode default>) -> PairAlignment

Optimal alignment of two nucleotide sequences (strings or encoded vectors)
with affine gap penalties, without band. `mode = :global` aligns both
sequences end to end (with the end gap penalties of `scoring`); `mode =
:query` fits the whole query into the target: gaps at the ends of the query
(unaligned target bases before and after) are free, so `target_start` and
`target_end` are the coordinates of the query on the target. See
[`Scoring`](@ref) for the defaults.
"""
function align(query, target; mode::Symbol = :global,
               scoring::Scoring = _default_scoring(_check_mode(mode)))
    _check_mode(mode)
    a = encode(query); b = encode(target)
    la = length(a); lb = length(b)
    (la == 0 || lb == 0) && throw(ArgumentError("sequences must not be empty"))
    ws = _Workspace()
    st = _align(ws, a, b, -la, lb, _kernel_scoring(scoring, mode))
    PairAlignment(st.score, st.matches, st.mismatches, st.internal_gap_openings, st.columns,
                  st.astart, st.aend, st.bstart, st.bend, _cigar(ws.ops, st, la))
end
