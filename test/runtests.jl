using NucleotideSearch
using NucleotideSearch: _Scoring, _Workspace, _align, _kernel_scoring
using Random
using Test

# Reference affine-gap DP with full matrices and the same cost model as the
# kernel: horizontal gaps (consuming `b`) in the first/last row use the end
# gap penalties of `a`, vertical gaps in the first/last column those of `b`.
function reference_score(a::Vector{UInt8}, b::Vector{UInt8}, sc::_Scoring)
    la = length(a); lb = length(b)
    neg = typemin(Int) ÷ 4
    H = fill(neg, la + 1, lb + 1); E = fill(neg, la + 1, lb + 1); F = fill(neg, la + 1, lb + 1)
    H[1, 1] = 0
    for i in 0:la, j in 0:lb
        (i == 0 && j == 0) && continue
        if j > 0
            row = i == 0 || i == la
            o = row ? sc.a_end_open : sc.gap_open; e = row ? sc.a_end_extend : sc.gap_extend
            E[i+1, j+1] = max(H[i+1, j] - o - e, E[i+1, j] - e)
        end
        if i > 0
            col = j == 0 || j == lb
            o = col ? sc.b_end_open : sc.gap_open; e = col ? sc.b_end_extend : sc.gap_extend
            F[i+1, j+1] = max(H[i, j+1] - o - e, F[i, j+1] - e)
        end
        diag = i > 0 && j > 0 ?
               H[i, j] + (a[i] == b[j] && a[i] < 0x04 ? sc.match : sc.mismatch) : neg
        H[i+1, j+1] = max(diag, E[i+1, j+1], F[i+1, j+1])
    end
    H[la+1, lb+1]
end

# Score of the path left in the workspace by `_align` (operations in reverse).
function path_score(ops, a, b, sc::_Scoring)
    la, lb = length(a), length(b)
    i = j = 0
    s = 0
    prev = 0x00
    for o in Iterators.reverse(ops)
        if o == UInt8('M')
            i += 1; j += 1
            s += (a[i] == b[j] && a[i] < 4) ? sc.match : sc.mismatch
        elseif o == UInt8('D')
            term = i == 0 || i == la
            s -= (term ? sc.a_end_extend : sc.gap_extend) +
                 (prev == o ? 0 : (term ? sc.a_end_open : sc.gap_open))
            j += 1
        else
            term = j == 0 || j == lb
            s -= (term ? sc.b_end_extend : sc.gap_extend) +
                 (prev == o ? 0 : (term ? sc.b_end_open : sc.gap_open))
            i += 1
        end
        prev = o
    end
    @assert i == la && j == lb
    s
end

randseq(rng, n) = rand(rng, UInt8(0):UInt8(3), n)

function mutate(rng, s::Vector{UInt8}, rate)
    out = UInt8[]
    for x in s
        r = rand(rng)
        if r < rate / 3
            push!(out, (x + rand(rng, 1:3)) % 4)       # substitution
        elseif r < 2rate / 3
            continue                                  # deletion
        elseif r < rate
            push!(out, x, rand(rng, UInt8(0):UInt8(3)))  # insertion
        else
            push!(out, x)
        end
    end
    out
end

decode(s) = String([b"ACGTN"[x+1] for x in s])

# Columns of a CIGAR string: (query bases, target bases).
function cigar_lengths(c::String)
    q = t = 0
    for m in eachmatch(r"(\d+)([MIDS])", c)
        n = parse(Int, m[1]); op = m[2][1]
        op in ('M', 'I', 'S') && (q += n)
        op in ('M', 'D') && (t += n)
    end
    q, t
end

@testset "NucleotideSearch" begin

@testset "encoding" begin
    @test encode("ACGTUacgtuN-R") == UInt8[0, 1, 2, 3, 3, 0, 1, 2, 3, 3, 4, 4, 4]
    @test reverse_complement(encode("AACGN")) == encode("NCGTT")
end

@testset "kernel against reference DP" begin
    rng = MersenneTwister(1)
    for trial in 1:400
        la = rand(rng, 1:60); lb = rand(rng, 1:60)
        a = randseq(rng, la)
        b = rand(rng) < 0.5 ? mutate(rng, a, 0.2) : randseq(rng, lb)
        isempty(b) && (b = randseq(rng, 3))
        rand(rng) < 0.1 && (a[rand(rng, 1:length(a))] = 0x04)
        sc = _Scoring(rand(rng, 1:3), -rand(rng, 1:5), rand(rng, 0:20), rand(rng, 0:3),
                      rand(rng, 0:5), rand(rng, 0:2), rand(rng, 0:5), rand(rng, 0:2))
        ws = _Workspace()
        st = _align(ws, a, b, -length(a), length(b), sc, Int16)
        @test st.score == reference_score(a, b, sc)
        @test path_score(ws.ops, a, b, sc) == st.score
        @test st.columns == length(ws.ops)
        @test st.matches + st.mismatches <= min(length(a), length(b))
        @test !st.touches_band_edge
        # 16- and 32-bit scores give the same alignment.
        @test _align(ws, a, b, -length(a), length(b), sc, Int32) == st
    end
end

@testset "kernel details" begin
    rng = MersenneTwister(42)
    sc = _Scoring(2, -4, 20, 2, 2, 1, 2, 1)
    ws = _Workspace()
    @test NucleotideSearch._score_type(300, 300, sc) === Int16
    @test NucleotideSearch._score_type(5_000, 5_000, sc) === Int32
    # Long sequences (32-bit scores).
    long = randseq(rng, 3_000)
    lb = mutate(rng, long, 0.01)
    st32 = _align(ws, long, lb, -length(long), length(lb), sc)
    @test st32.score == path_score(ws.ops, long, lb, sc)
    @test st32.matches >= 2_900
    a = encode("ACGTACGTAC"); b = encode("ACGTTCGTAC")
    st = _align(ws, a, b, -2, 2, sc)
    @test (st.matches, st.mismatches, st.columns, st.terminal_gap_columns, st.gap_openings) ==
          (9, 1, 10, 0, 0)
    @test NucleotideSearch._identity(st, 10, 10, 1) == 90
    # End gaps: b is a with 5 extra bases at the end.
    a = encode("ACGTACGTACGTACGTACGT"); b = encode("ACGTACGTACGTACGTACGTGGGGG")
    st = _align(ws, a, b, -20, 25, sc)
    @test (st.matches, st.columns, st.terminal_gap_columns, st.gap_openings) == (20, 25, 5, 1)
    @test (st.astart, st.aend, st.bstart, st.bend) == (1, 20, 1, 20)
    @test NucleotideSearch._identity(st, 20, 25, 0) == 100
    @test NucleotideSearch._identity(st, 20, 25, 1) == 80
    @test NucleotideSearch._identity(st, 20, 25, 2) == 100
    @test NucleotideSearch._identity(st, 20, 25, 3) ≈ 96
    @test NucleotideSearch._identity(st, 20, 25, 4) == 80
    # N never matches.
    st = _align(ws, encode("ANGT"), encode("ANGT"), -4, 4, sc)
    @test st.matches == 3 && st.mismatches == 1
end

@testset "align: global and query modes" begin
    rng = MersenneTwister(2)
    target = randseq(rng, 300)
    q = target[101:200]
    r = align(q, target; mode = :query)
    @test (r.target_start, r.target_end) == (101, 200)
    @test (r.query_start, r.query_end) == (1, 100)
    @test r.matches == 100 && r.mismatches == 0
    @test r.cigar == "100M"

    # Mismatches at both ends of the query are not clipped.
    q2 = copy(q); q2[1] = (q2[1] + 1) % 4; q2[2] = (q2[2] + 1) % 4; q2[end] = (q2[end] + 2) % 4
    r2 = align(q2, target; mode = :query)
    @test (r2.target_start, r2.target_end) == (101, 200)
    @test (r2.query_start, r2.query_end) == (1, 100)
    @test r2.mismatches == 3

    # An indel inside the query keeps the coordinates exact. (Within the first
    # few bases, shifting the start can score better than opening a gap.)
    q3 = vcat(q[1:29], q[31:end])
    r3 = align(q3, target; mode = :query)
    @test (r3.target_start, r3.target_end) == (101, 200)
    @test cigar_lengths(r3.cigar) == (99, 100)
    @test r3.score == reference_score(encode(decode(q3)), target,
                                      _kernel_scoring(NucleotideSearch._default_scoring(:query), :query))

    # Query longer than a truncated target: unaligned query bases are soft clips.
    r4 = align(q, target[1:150]; mode = :query)
    @test (r4.target_start, r4.target_end) == (101, 150)
    @test (r4.query_start, r4.query_end) == (1, 50)
    @test r4.cigar == "50M50S"

    g = align("ACGTACGT", "ACGTTCGT")
    @test g.matches == 7 && g.mismatches == 1 && g.cigar == "8M"
    @test_throws ArgumentError align("", "ACGT")
    @test_throws ArgumentError align("ACGT", "ACGT"; mode = :local)
    @test_throws ArgumentError Scoring(gap_open = -1)
end

@testset "database: build, save, load" begin
    db = SequenceDB(["a", "b", "c"], ["ACGTACGTAC", "TTTTGGGGCCCC", "ACGTNACGTA"]; k = 4)
    @test nsequences(db) == 3
    path = tempname()
    try
        save_db(path, db)
        db2 = load_db(path)
        for f in fieldnames(SequenceDB)
            @test getfield(db, f) == getfield(db2, f)
        end
        write(path, "nothing")
        @test_throws ArgumentError load_db(path)
    finally
        rm(path; force = true)
    end
    fasta = tempname()
    try
        write(fasta, ">a desc\nACGTAC\nGTAC\n\n>b\nTTTTGGGGCCCC\n>c\nACGTNACGTA\n")
        @test SequenceDB(fasta; k = 4).data == db.data
    finally
        rm(fasta; force = true)
    end
    @test_throws ArgumentError SequenceDB(["a"], ["ACGT"]; k = 2)
    @test_throws DimensionMismatch SequenceDB(["a"], ["ACGT", "ACGT"])
end

@testset "search: coordinates on a reference database" begin
    rng = MersenneTwister(3)
    refs = [randseq(rng, rand(rng, 600:900)) for _ in 1:200]
    db = SequenceDB(["ref$i" for i in eachindex(refs)], refs)
    queries = Vector{UInt8}[]; truth = Tuple{Int,Char,Int,Int}[]
    for i in 1:60
        t = rand(rng, eachindex(refs)); start = rand(rng, 1:300); len = rand(rng, 150:300)
        stop = start + len - 1
        q = copy(refs[t][start:stop])
        for p in rand(rng, 1:len, 8)       # substitutions, including near the ends
            q[p] = (q[p] + 1) % 4
        end
        q[1] = (q[1] + 2) % 4; q[end] = (q[end] + 2) % 4
        strand = rand(rng) < 0.5 ? '+' : '-'
        push!(queries, strand == '+' ? q : reverse_complement(q))
        push!(truth, (t, strand, start, stop))
    end
    hits = search(queries, db; cigar = true)
    for (h, (t, strand, start, stop), q) in zip(hits, truth, queries)
        @test !isempty(h)
        best = h[1]
        @test best.target == t
        @test best.strand == strand
        @test (best.target_start, best.target_end) == (start, stop)
        @test (best.query_start, best.query_end) == (1, length(q))
        @test cigar_lengths(best.cigar) == (length(q), stop - start + 1)
        # Same result as the unbanded alignment on the whole target.
        a = strand == '+' ? q : reverse_complement(q)
        full = align(a, refs[t]; mode = :query)
        @test best.score == full.score
    end
    @test search(queries, db; strand = :plus)[findfirst(x -> x[2] == '-', truth)] |>
          h -> all(x -> x.strand == '+', h)
    # One query given as a string.
    @test search(decode(queries[1]), db; cigar = true) == hits[1]
end

@testset "search: indels and window fallback" begin
    rng = MersenneTwister(4)
    refs = [randseq(rng, 1500) for _ in 1:20]
    db = SequenceDB(["r$i" for i in 1:20], refs)
    for trial in 1:40
        t = rand(rng, 1:20); start = rand(rng, 1:1200); stop = start + 249
        q = mutate(rng, refs[t][start:stop], 0.06)
        h = search(q, db; strand = :plus, min_identity = 70)
        isempty(h) && continue
        full = align(q, refs[t]; mode = :query)
        best = only(filter(x -> x.target == t, h))
        @test best.score == full.score
        @test (best.target_start, best.target_end) == (full.target_start, full.target_end)
    end
end

@testset "search_self" begin
    rng = MersenneTwister(5)
    base = randseq(rng, 300)
    seqs = [[base]; [mutate(rng, base, 0.03) for _ in 1:15]; [randseq(rng, 300) for _ in 1:10]]
    db = SequenceDB(["s$i" for i in eachindex(seqs)], seqs)
    res = search_self(db; maxhits = 0)
    @test all(isempty, res[17:end])
    @test all(h -> all(x -> x.target <= 16 && x.target != x.query, h), res[1:16])
    # Symmetric identities: (q, t) and (t, q) share one alignment.
    ids = Dict((x.query, x.target) => x.identity for h in res for x in h)
    @test all(((q, t),) -> get(ids, (t, q), ids[(q, t)]) == ids[(q, t)], keys(ids))
    # Coordinates refer to the query of each hit.
    for h in res, x in h
        @test x.query_end <= length(seqs[x.query]) && x.target_end <= length(seqs[x.target])
    end
    exhaustive = search_self(db; maxrejects = 0, minwordmatches = 0, maxhits = 0)
    @test exhaustive == res
    @test_throws ArgumentError search_self(db; iddef = 5)
    @test_throws ArgumentError search_self(db; query_cov = 2)
    @test_logs (:warn, r"percentage") search_self(db; min_identity = 0.9)
end

end
