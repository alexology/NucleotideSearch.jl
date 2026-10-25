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
    max_postings::Int       # k-mers in more targets are not counted (`search` only)
end

function _settings(mode, min_identity, iddef, query_cov, maxhits, maxaccepts, maxrejects,
                   minwordmatches, scoring::Scoring, cigar; max_postings::Int = typemax(Int))
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
                    _kernel_scoring(scoring, mode), cigar, max_postings)
end

# Diagonal band containing every alignment of sequences of lengths `la` and
# `lb` that can be accepted. Under iddef 1/4 an alignment with identity >= t
# has at most `slack` gap columns, so it stays within (slack -/+ (lb - la)) / 2
# of the main diagonal. The other definitions ignore end gaps (0, 2) or count
# only their openings (3), so the overhangs are bounded by the number `A` of
# aligned pairs instead: every point of the path has between `la - A` unpaired
# bases of `a` and `lb - A` of `b` behind it, so `A - la <= j - i <= lb - A`.
# An accepted hit has `A >= query_cov` times the query length (either sequence
# in `search_self`, hence the shorter one) and, under iddef 0, `A >= matches >=
# t * min(la, lb)`.
function _band(la::Int, lb::Int, s::_SearchSettings)
    t = s.min_identity
    delta = lb - la
    if s.iddef == 1 || s.iddef == 4
        slack = t > 0 ? floor(Int, (1 - t) / t * min(la, lb)) : la + lb
        max(-la, min(0, delta, -fld(slack - delta, 2))), min(lb, max(0, delta, fld(slack + delta, 2)))
    else
        short = min(la, lb)
        amin = s.query_cov * short
        s.iddef == 0 && (amin = max(amin, s.id_threshold / 100 * short))
        # One pair less, for rounding
        amin = clamp(floor(Int, amin) - 1, 0, short)
        max(-la, amin - la), min(lb, lb - amin)
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
    qkp::NTuple{2,Vector{UInt32}}  # per strand: k-mers occurring once in the query,
    qpos::NTuple{2,Vector{Int32}}  # sorted, with their end positions
    diags::Vector{Int}
end

_SearchState(n::Int) = _SearchState(_Workspace(), zeros(Int32, n), Int[], Int[],
                                    (UInt32[], UInt32[]), (Int32[], Int32[]), Int[])

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
"""
    SearchReport

How the search of one query ended, from [`search`](@ref) and
[`search_self`](@ref) with `reports`. Fields:

- `candidates`: targets (× strands) sharing enough k-mers with the query to
  be considered (`minwordmatches`);
- `examined`: candidates checked, in ranking order (length test, then
  alignment);
- `accepted`, `rejected`: of the examined candidates;
- `stop`: `:maxaccepts` or `:maxrejects` if that limit ended the search,
  even on the last candidate; `:exhausted` if candidates ran out without
  reaching either limit;
- `exhaustive`: every target of the database (every strand) was examined, so
  no better hit can exist elsewhere, independently of `stop` (the query
  itself is excluded in `search_self`). Otherwise the hits are the best among
  the examined candidates only; each hit's alignment is optimal for its pair
  either way.
"""
struct SearchReport
    candidates::Int
    examined::Int
    accepted::Int
    rejected::Int
    stop::Symbol
    exhaustive::Bool
end

function _report(candidates, examined, accepted, rejected, stop, possible)
    SearchReport(candidates, examined, accepted, rejected, stop,
                 examined == possible)
end

_check_reports(reports, n) = reports === nothing || resize!(reports, n)

function _finish!(hits::Vector{Hit}, s::_SearchSettings)
    sort!(hits; alg = MergeSort, by = h -> -h.identity)
    s.maxhits > 0 && length(hits) > s.maxhits && resize!(hits, s.maxhits)
    hits
end

# `f(state, q)` for every query `q` on at most `threads` tasks, each with its
# own search state (scratch sized by the database); `threads = 1` runs on the
# calling task. Queries are handed out in chunks, so the result does not
# depend on the scheduling.
function _run_parallel(f, n::Int, ntargets::Int, threads::Int)
    threads >= 1 || throw(ArgumentError("threads must be at least 1"))
    nworkers = min(threads, n)
    if nworkers <= 1
        state = _SearchState(ntargets)
        for q in 1:n
            f(state, q)
        end
        return
    end
    chunk = max(1, cld(n, 8 * nworkers))
    next = Threads.Atomic{Int}(1)
    tasks = map(1:nworkers) do _
        Threads.@spawn begin
            # `local`: the serial branch's `state` must not be captured and shared
            local state = _SearchState(ntargets)
            while (lo = Threads.atomic_add!(next, chunk)) <= n
                for q in lo:min(n, lo + chunk - 1)
                    f(state, q)
                end
            end
        end
    end
    failure = nothing
    for task in tasks
        try
            wait(task)
        catch e
            failure === nothing && (failure = e)
        end
    end
    failure === nothing || throw(failure)
end

# ---------------------------------------------------------------------------
# All against all

# Alignments shared between the two directions of a pair, in lock-striped
# dictionaries. An entry is removed when it is used the second time. Each
# stripe has a fixed share of the entry budget; full stripes skip insertion.
struct _PairCache
    tables::Vector{Dict{Int,_AlignStats}}
    locks::Vector{ReentrantLock}
    max_entries::Int
end

function _PairCache(nstripes::Int; max_entries::Int = 100_000)
    nstripes > 0 || throw(ArgumentError("cache stripe count must be positive"))
    max_entries >= 0 || throw(ArgumentError("cache_entries must be non-negative"))
    _PairCache([Dict{Int,_AlignStats}() for _ in 1:nstripes],
               [ReentrantLock() for _ in 1:nstripes], max_entries)
end

function _get_or_align!(compute, cache::_PairCache, key::Int)
    stripe = mod1(key, length(cache.tables))
    table = cache.tables[stripe]
    lk = cache.locks[stripe]
    cached = lock(() -> pop!(table, key, nothing), lk)
    cached === nothing || return cached
    st = compute()   # outside the lock; a concurrent duplicate gives the same result
    lock(lk) do
        limit = div(cache.max_entries, length(cache.tables)) +
                Int(stripe <= rem(cache.max_entries, length(cache.tables)))
        if haskey(table, key)
            delete!(table, key)
        elseif length(table) < limit
            table[key] = st
        end
    end
    st
end

# With the same end gap penalties for both sequences, swapping them gives the
# same optimal alignment.
_symmetric(sc::_Scoring) = sc.a_end_open == sc.b_end_open && sc.a_end_extend == sc.b_end_extend

# Alignment of sequences `x` and `y`. With symmetric scoring a pair is always
# aligned in the canonical orientation `x < y`, so (q, t) and (t, q) share one
# alignment (identity and number of aligned pairs are symmetric); otherwise each
# direction is aligned on its own.
function _pair_stats(ws::_Workspace, db::SequenceDB, x::Int, y::Int, s::_SearchSettings)
    a = _seq(db, x); b = _seq(db, y)
    _global_stats(ws, a, b, s)
end

function _global_stats(ws::_Workspace, a, b, s::_SearchSettings)
    la = length(a); lb = length(b)
    klo, khi = _band(la, lb, s)
    st = _align(ws, a, b, klo, khi, s.scoring)
    # A tie outside the band can change identity or coverage under the kernel's
    # traceback rules, so only a strict score bound certifies the result.
    if ws.outside_bound >= st.score
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
    candidates = _rank_candidates!(st, _unique_kmers(_seq(db, q), db.k), db, s)
    ncand = length(candidates) - (q in candidates)
    stop = :exhausted
    for tg in candidates
        tg == q && continue
        lb = _seqlen(db, tg)
        ok = false
        if _length_compatible(la, lb, s)
            if _symmetric(s.scoring)
                x, y = minmax(q, tg)
                stats = _get_or_align!(cache, (x - 1) * n + y) do
                    _pair_stats(st.ws, db, x, y, s)
                end
                q > tg && (stats = _swap(stats))
            else
                stats = _pair_stats(st.ws, db, q, tg, s)
            end
            id = _identity(stats, la, lb, s.iddef)
            ok = id >= s.id_threshold && stats.matches + stats.mismatches >= s.query_cov * la
            ok && push!(accepted, _hit(q, tg, '+', id, stats, la, ""))
        end
        if ok
            s.maxaccepts > 0 && length(accepted) >= s.maxaccepts && (stop = :maxaccepts; break)
        else
            rejected += 1
            s.maxrejects > 0 && rejected >= s.maxrejects && (stop = :maxrejects; break)
        end
    end
    report = _report(ncand, length(accepted) + rejected, length(accepted), rejected, stop, n - 1)
    _finish!(accepted, s), report
end

"""
    search_self(db::SequenceDB; min_identity = 84, iddef = 1, query_cov = 0.9,
                maxhits = 10, maxaccepts = 0, maxrejects = 32,
                minwordmatches = 12, scoring = Scoring(), reports = nothing,
                threads = Threads.nthreads(), cache_entries = 100_000)
        -> Vector{Vector{Hit}}

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
(exhaustive search; slow for large data sets). Pass `reports =
SearchReport[]` to receive, for each query, a [`SearchReport`](@ref): how
many candidates were examined, why the search stopped and whether it covered
the whole database.

Identity definitions (`iddef`, as in VSEARCH): `0`: matches / length of the
shorter sequence; `1` and `4`: matches / alignment columns; `2`: matches /
alignment columns excluding end gaps; `3`: `1 - (mismatches + gap openings) /
length of the longer sequence`. Bases other than A, C, G, T/U never count as
matches.

The alignment is restricted to a diagonal band wide enough for every
alignment that can be accepted (identity under `iddef` and `query_cov`). A
score bound certifies that no path outside the band can tie or beat its
optimum; otherwise the pair is recomputed without band. The reciprocal cache
holds at most `cache_entries` alignments (0 disables it); full stripes
recompute pairs instead, preserving the same results. Queries run in
parallel on `threads` tasks (default `Threads.nthreads()`; `1` runs on the
calling task, e.g. inside an already parallel pipeline); the result does not
depend on the number of threads.
"""
function search_self(db::SequenceDB; min_identity::Real = 84, iddef::Integer = 1,
                     query_cov::Real = 0.9, maxhits::Integer = 10, maxaccepts::Integer = 0,
                     maxrejects::Integer = 32, minwordmatches::Integer = 12,
                     scoring::Scoring = Scoring(),
                     reports::Union{Nothing,Vector{SearchReport}} = nothing,
                     threads::Integer = Threads.nthreads(), cache_entries::Integer = 100_000)
    s = _settings(:global, min_identity, iddef, query_cov, maxhits, maxaccepts, maxrejects,
                  minwordmatches, scoring, false)
    n = nsequences(db)
    cache = _PairCache(64; max_entries = Int(cache_entries))
    results = Vector{Vector{Hit}}(undef, n)
    _check_reports(reports, n)
    _run_parallel(n, n, Int(threads)) do state, q
        results[q], report = _search_self_query(state, q, db, cache, s)
        reports === nothing || (reports[q] = report)
    end
    results
end

# ---------------------------------------------------------------------------
# Queries against a database

# Most frequent diagonal (target position - query position) among the k-mers
# that occur once in the query and appear in the target; `nothing` if none.
function _seed_diagonal(st::_SearchState, si::Int, t::AbstractVector{UInt8}, k::Int)
    kp = st.qkp[si]; pos = st.qpos[si]
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

# Proof that an alignment of the query (length `la`) scoring `score` in a
# window of the target is optimal on the whole target (q-gram lemma). An
# alignment scoring at least as much loses at most `deficit` points against a perfect
# match of the query, and each difference destroys a bounded number of query
# k-mers per point lost: a mismatch at most k, a run of n inserted query bases
# n + k - 1, a run of deleted target bases k - 1, an unaligned query end of n
# bases n. Such an alignment therefore keeps at least `need` of the `nunique`
# k-mers unique in the query as seeds (`diags`, sorted), and its diagonals span
# at most `width` (its inserted and deleted bases). Every group of `need` seeds
# that close must lie where the window alignment could have used them: inside
# the window band (`kwlo:kwhi`, whole-target diagonals) and the window `lo:hi`.
function _window_optimal(diags::Vector{Int}, nunique::Int, score::Int, la::Int, lb::Int,
                         lo::Int, hi::Int, kwlo::Int, kwhi::Int, k::Int, sc::_Scoring)
    m = Int(sc.match); x = m - Int(sc.mismatch)
    go = Int(sc.gap_open); ge = Int(sc.gap_extend)
    ee = Int(sc.b_end_extend)   # unaligned query ends
    min(m, x, ge, m + ee) > 0 || return false
    rate = max(k / x, k / (m + go + ge), 1 / (m + ge), (k - 1) / (go + ge), 1 / (m + ee))
    # Ties can change identity, coverage and coordinates under kernel tie rules.
    deficit = m * la - score
    deficit < 0 && return true
    need = nunique - floor(Int, rate * deficit)
    need > 0 || return false
    width = fld(deficit, ge)
    n = length(diags); j = 1
    for i in 1:n
        j = max(j, i)
        while j < n && diags[j+1] - diags[i] <= width
            j += 1
        end
        j - i + 1 >= need || continue
        # `d` is the lowest seed diagonal of a possible equal or better alignment
        d = diags[i]
        kwlo <= d - width && d + width <= kwhi || return false
        lo == 1 || d - width + 1 >= lo || return false
        hi == lb || d + width + la <= hi || return false
    end
    true
end

# Fitting (`:query` mode) of query `a` into target `b`, in a window of the
# target around the seed diagonal. The result is recomputed on the whole target
# without band if the path touches the band edge or a window boundary inside
# the target, or if an equal or better alignment elsewhere cannot be ruled out
# (`_window_optimal`). Coordinates refer to the whole target.
function _query_stats(st::_SearchState, si::Int, a::AbstractVector{UInt8}, b::AbstractVector{UInt8},
                      k::Int, s::_SearchSettings)
    la = length(a); lb = length(b)
    t = s.min_identity
    slack = t > 0 ? min(la + lb, floor(Int, (1 - t) / t * la)) : la + lb
    d = _seed_diagonal(st, si, b, k)
    if d !== nothing
        margin = slack + 16
        lo = max(1, d + 1 - margin); hi = min(lb, d + la + margin)
        if lo <= hi
            w = view(b, lo:hi); lw = length(w)
            klo = max(-la, -slack); khi = min(lw, lw - la + slack)
            if klo <= min(0, lw - la) && khi >= max(0, lw - la)
                r = _align(st.ws, a, w, klo, khi, s.scoring)
                inside = r.aend > 0 && (lo == 1 || r.bstart > 1) && (hi == lb || r.bend < lw)
                if inside && !r.touches_band_edge &&
                   _window_optimal(st.diags, length(st.qkp[si]), r.score, la, lb, lo, hi,
                                   klo + lo - 1, khi + lo - 1, k, s.scoring)
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

function _target_stats(st::_SearchState, si::Int, a, b, k::Int, s::_SearchSettings)
    if s.mode === :query
        _query_stats(st, si, a, b, k, s)
    else
        _global_stats(st.ws, a, b, s)
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
        s.mode === :query && _kmer_positions!(st.qkp[si], st.qpos[si], a, k)
    end
    accepted = Hit[]
    rejected = 0
    # The candidates are ordered in blocks (partial sort): usually only the
    # first few hundred are aligned before `maxaccepts`/`maxrejects` stop the
    # search, out of up to all references. The order is the full sort order.
    n = length(cands); lo = 1; done = false; stop = :exhausted
    while lo <= n && !done
        hi = min(n, lo + 255)
        partialsort!(view(cands, lo:n), 1:(hi - lo + 1))
        for c in lo:hi
            _, tg, si = cands[c]
            b = _seq(db, tg)
            lb = length(b)
            a = oriented[si]
            ok = false
            if _length_compatible(la, lb, s)
                r = _target_stats(st, si, a, b, k, s)
                id = _identity(r, la, lb, s.iddef)
                ok = r.aend > 0 && id >= s.id_threshold &&
                     r.matches + r.mismatches >= s.query_cov * la
                # The workspace still holds the path of `r`.
                ok && push!(accepted, _hit(q, tg, strands[si], id, r, la,
                                           s.cigar ? _cigar(st.ws.ops, r, la) : ""))
            end
            if ok
                s.maxaccepts > 0 && length(accepted) >= s.maxaccepts && (done = true; stop = :maxaccepts; break)
            else
                rejected += 1
                s.maxrejects > 0 && rejected >= s.maxrejects && (done = true; stop = :maxrejects; break)
            end
        end
        lo = hi + 1
    end
    report = _report(n, length(accepted) + rejected, length(accepted), rejected, stop,
                     nsequences(db) * length(strands))
    _finish!(accepted, s), report
end

# K-mers in more than this many targets are not counted: a fraction of the
# database, but never fewer than 50,000 targets. On databases of 5,000-10,000
# references a cap made the best hits worse (fewer good candidates) for little
# speed; on SILVA (313,502 references) a cap of 0.2 halved the search time
# with the same results.
function _max_postings(frac::Real, n::Int)
    0 < frac <= 1 || throw(ArgumentError("max_kmer_frequency must be in (0, 1]"))
    frac == 1 ? typemax(Int) : max(50_000, floor(Int, frac * n))
end

# Appends (-shared k-mers, target, `si`) for the targets sharing enough of the
# query k-mers `qk` (all targets if the threshold is 0). K-mers found in more
# than `s.max_postings` targets are not counted, unless fewer than
# `minwordmatches` query k-mers would remain.
function _collect_candidates!(cands, st::_SearchState, qk::Vector{UInt32}, db::SequenceDB,
                              s::_SearchSettings, si::Int)
    shared = st.shared; touched = empty!(st.touched)
    offsets = db.offsets; postings = db.postings
    if s.max_postings < length(postings)
        rare = filter(w -> offsets[w+2] - offsets[w+1] <= s.max_postings, qk)
        length(rare) >= s.minwordmatches && (qk = rare)
    end
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
           scoring = <mode default>, cigar = false, max_kmer_frequency = 1,
           reports = nothing, threads = Threads.nthreads())
        -> Vector{Vector{Hit}}

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
  whole target if it reaches the window or band boundaries, or if the shared
  k-mers do not rule out a better alignment elsewhere in the target.
- `mode = :global`: end-to-end alignment as in [`search_self`](@ref).
- `strand = :both` also searches the reverse complement of each query
  (`:plus`: forward only).

Candidate selection, acceptance (`min_identity`, `query_cov`), `maxaccepts`,
`maxrejects`, `maxhits` and `iddef` are as in [`search_self`](@ref); `iddef =
2` (default here) ignores the unaligned ends of the target. Candidates of both
strands are ranked together. `cigar = true` adds the CIGAR string of each hit.
`reports = SearchReport[]` is filled with one [`SearchReport`](@ref) per query.
Queries run on `threads` tasks, as in [`search_self`](@ref).

`max_kmer_frequency < 1` skips, when counting shared k-mers, the k-mers found
in more than that fraction of the references and in more than 50,000 of
them (conserved regions of large databases, which make every reference a
candidate and dominate the search time); if fewer than `minwordmatches`
query k-mers remain, all are used. On SILVA (313,502 16S references)
`max_kmer_frequency = 0.2` halved the search time with the same hits.
"""
function search(queries::AbstractVector, db::SequenceDB; mode::Symbol = :query,
                strand::Symbol = :both, min_identity::Real = 80, iddef::Integer = 2,
                query_cov::Real = 0.9, maxhits::Integer = 10, maxaccepts::Integer = 10,
                maxrejects::Integer = 32, minwordmatches::Integer = 12,
                scoring::Scoring = _default_scoring(_check_mode(mode)), cigar::Bool = false,
                max_kmer_frequency::Real = 1.0,
                reports::Union{Nothing,Vector{SearchReport}} = nothing,
                threads::Integer = Threads.nthreads())
    _check_mode(mode)
    strand in (:both, :plus) || throw(ArgumentError("strand must be :both or :plus"))
    # Columns and gap openings that include end gaps depend on the target window
    # in `:query` mode, so only iddef 0 and 2 are meaningful there.
    mode === :query && !(iddef in (0, 2)) &&
        throw(ArgumentError("in :query mode iddef must be 0 or 2"))
    strands = strand === :both ? ('+', '-') : ('+',)
    s = _settings(mode, min_identity, iddef, query_cov, maxhits, maxaccepts, maxrejects,
                  minwordmatches, scoring, cigar;
                  max_postings = _max_postings(max_kmer_frequency, nsequences(db)))
    encoded = [collect(encode(x)) for x in queries]
    n = length(encoded)
    results = Vector{Vector{Hit}}(undef, n)
    _check_reports(reports, n)
    _run_parallel(n, nsequences(db), Int(threads)) do state, q
        results[q], report = _search_query(state, q, encoded[q], db, strands, s)
        reports === nothing || (reports[q] = report)
    end
    results
end

search(query::Union{AbstractString,AbstractVector{UInt8}}, db::SequenceDB; kwargs...) =
    only(search([query], db; kwargs...))
