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
        # Soft clipping: start anywhere after clipping the first i bases of a.
        sc.clip_open >= 0 && i >= 1 &&
            (H[i+1, j+1] = max(H[i+1, j+1], -(sc.clip_open + sc.clip_extend * i)))
    end
    best = H[la+1, lb+1]
    if sc.clip_open >= 0                    # or end anywhere, clipping the rest of a
        for i in 0:la, j in 0:lb
            best = max(best, H[i+1, j+1] - (sc.clip_open + sc.clip_extend * (la - i)))
        end
    end
    best
end

# Score of a clipped path: the part from the first to the last aligned pair
# with internal penalties, plus the clipping costs of the unaligned ends of `a`.
function clipped_path_score(ops, a, b, sc::_Scoring, st)
    st.aend == 0 && return -(sc.clip_open + sc.clip_extend * length(a))
    i = st.astart - 1; j = st.bstart - 1
    s = 0; prev = 0x00
    lo = findfirst(==(UInt8('M')), ops); hi = findlast(==(UInt8('M')), ops)
    for k in hi:-1:lo
        o = ops[k]
        if o == UInt8('M')
            i += 1; j += 1
            s += (a[i] == b[j] && a[i] < 4) ? sc.match : sc.mismatch
        else
            s -= sc.gap_extend + (prev == o ? 0 : sc.gap_open)
            o == UInt8('D') ? (j += 1) : (i += 1)
        end
        prev = o
    end
    clip(n) = n == 0 ? 0 : sc.clip_open + sc.clip_extend * n
    s - clip(st.astart - 1) - clip(length(a) - st.aend)
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
        clip = rand(rng) < 0.4
        if clip
            # `:query` mode: clipping costs as much as target end gaps. Gaps
            # cost strictly more than clipping (open > clip open, extend >= clip
            # extend), as with the defaults (20 + 2n against 20 + n), so that an
            # unaligned end is always a clip and `clipped_path_score` applies.
            go = rand(rng, 1:20); co = rand(rng, 0:go-1); ce = rand(rng, 0:2)
            sc = _Scoring(rand(rng, 1:3), -rand(rng, 1:5), go, rand(rng, ce:3),
                          0, 0, co, ce, co, ce)
        else
            sc = _Scoring(rand(rng, 1:3), -rand(rng, 1:5), rand(rng, 0:20), rand(rng, 0:3),
                          rand(rng, 0:5), rand(rng, 0:2), rand(rng, 0:5), rand(rng, 0:2))
        end
        ws = _Workspace()
        st = _align(ws, a, b, -length(a), length(b), sc, Int16)
        @test st.score == reference_score(a, b, sc)
        if clip
            @test clipped_path_score(ws.ops, a, b, sc, st) == st.score
        else
            @test path_score(ws.ops, a, b, sc) == st.score
        end
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
    @test r.gap_openings == 0          # the unaligned target ends are not gaps

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

    # An end that differs entirely while the target continues is soft-clipped
    # (20 + n) instead of being forced into the alignment.
    q5 = copy(q); q5[71:100] .= (q5[71:100] .+ 0x01) .% 0x04
    r5 = align(q5, target; mode = :query)
    @test (r5.query_start, r5.query_end) == (1, 70)
    @test (r5.target_start, r5.target_end) == (101, 170)
    @test r5.cigar == "70M30S"
    q6 = copy(q); q6[1:25] .= (q6[1:25] .+ 0x02) .% 0x04
    r6 = align(q6, target; mode = :query)
    @test r6.cigar == "25S75M" && (r6.target_start, r6.target_end) == (126, 200)
    # Global mode never clips.
    @test !NucleotideSearch._clipping(_kernel_scoring(Scoring(), :global))
    @test NucleotideSearch._clipping(_kernel_scoring(Scoring(), :query))

    g = align("ACGTACGT", "ACGTTCGT")
    @test g.matches == 7 && g.mismatches == 1 && g.cigar == "8M"
    @test_throws ArgumentError align("", "ACGT")
    @test_throws ArgumentError align("ACGT", "ACGT"; mode = :local)
    @test_throws ArgumentError Scoring(gap_open = -1)
end

@testset "database: reject corrupt storage before search" begin
    mktempdir() do dir
        path = joinpath(dir, "index.nsdb")
        valid = SequenceDB(["a", "b", "c"], ["AAAAAA", "AAAAAA", "CCCCCC"]; k = 3)
        for damage in (db -> (db.starts[1] = 0),
                       db -> (db.starts[2] = db.starts[3] + 1),
                       db -> (db.starts[end] += 1),
                       db -> (db.offsets[1] = 1),
                       db -> (db.offsets[2] = -1),
                       db -> (db.offsets[end] += 1),
                       db -> (db.postings[1] = 0),
                       db -> (db.postings[1] = 4),
                       db -> (db.postings[1] = db.postings[2]),
                       db -> reverse!(view(db.postings, 1:2)),
                       db -> (db.data[1] = 0xff),
                       db -> pop!(db.starts),
                       db -> pop!(db.offsets))
            damaged = deepcopy(valid)
            damage(damaged)
            save_db(path, damaged)
            @test_throws ArgumentError load_db(path)
        end
        save_db(path, valid)
        bytes = read(path)
        for nbytes in (0, 4, 7, 15, 23, length(bytes) - 1)
            write(path, bytes[1:nbytes])
            @test_throws ArgumentError load_db(path)
        end
        for (position, value) in ((8, Int64(0)), (8, Int64(13)), (16, Int64(-1)),
                                  (16, typemax(Int64)), (24, Int64(-1)), (24, typemax(Int64)))
            write(path, bytes)
            open(path, "r+") do io
                seek(io, position)
                write(io, htol(value))
            end
            @test_throws ArgumentError load_db(path)
        end
        write(path, [bytes; 0x00])
        @test_throws ArgumentError load_db(path)
        for db in (SequenceDB(String[], String[]; k = 3),
                   SequenceDB(["a", "b", "c"], ["", "NNN", ""]; k = 3), valid)
            save_db(path, db)
            loaded = load_db(path)
            @test all(getfield(db, f) == getfield(loaded, f) for f in fieldnames(SequenceDB))
        end
    end
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

@testset "read_fasta" begin
    fasta = tempname()
    readfa(text) = (write(fasta, text); NucleotideSearch.read_fasta(fasta))
    try
        @test readfa(">a desc\nACGT\nAC\n>b\nTTTT\n") == (["a", "b"], ["ACGTAC", "TTTT"])
        @test readfa(">a desc\r\nACGT\r\nAC\r\n>b\r\nTTTT\r\n") == (["a", "b"], ["ACGTAC", "TTTT"])
        @test readfa("\n>a\nACGT\n\nAC\n\n>b\nTTTT") == (["a", "b"], ["ACGTAC", "TTTT"])
        @test readfa(">a\nacgtNRY-.\n") == (["a"], ["acgtNRY-."])
        @test readfa(">a  \n  AC GT \n\tAC\n") == (["a"], ["ACGTAC"])
        @test readfa("> a desc\nACGT\n") == (["a"], ["ACGT"])
        @test readfa(">\nACGT\n") == ([""], ["ACGT"])
        @test readfa(">a\n\n>b\nACGT\n") == (["a", "b"], ["", "ACGT"])
        @test readfa("") == (String[], String[])
        @test_throws ArgumentError readfa("ACGT\n>a\nAC\n")
        @test_throws ArgumentError readfa(">a\n>b\nACGT\n")
    finally
        rm(fasta; force = true)
    end
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
        @test best.aligned_columns ==
              sum(parse(Int, m[1]) for m in eachmatch(r"(\d+)[MID]", best.cigar))
        @test best.gap_openings == count(m -> m[1] != "M", eachmatch(r"\d+([MID])", best.cigar))
        # Same result as the unbanded alignment on the whole target.
        a = strand == '+' ? q : reverse_complement(q)
        full = align(a, refs[t]; mode = :query)
        @test best.score == full.score
    end
    # On small databases the k-mer frequency cap (>= 50,000 targets) changes nothing.
    @test search(queries, db; max_kmer_frequency = 0.2) == search(queries, db)
    @test_throws ArgumentError search(queries, db; max_kmer_frequency = 0)
    @test_throws ArgumentError search(queries, db; max_kmer_frequency = 1.5)
    @test_throws ArgumentError search(queries, db; iddef = 1)
    @test_throws ArgumentError search(queries, db; strand = :minus)
    @test length(search(queries[1:2], db; mode = :global, iddef = 1, min_identity = 10,
                        query_cov = 0.1)) == 2
    @test search(queries, db; strand = :plus)[findfirst(x -> x[2] == '-', truth)] |>
          h -> all(x -> x.strand == '+', h)
    # One query given as a string.
    @test search(decode(queries[1]), db; cigar = true) == hits[1]
end

@testset "search: frequent k-mers are not counted" begin
    # A 40-base block shared by all targets; target 3 also shares the query's
    # unique part. With the cap, only the unique part counts.
    # Fixed sequences avoid accidental shared k-mers when Julia's RNG changes.
    block = encode("TTGAGCTCACTTGCCCGTCTCTAAGCATACCCCCGAGGAT")
    uniq = encode("TCTCTTCGCTGGCAAATAGCCGGCCACGCATGAATTGCAATGTAGCGCCGCTCGCTTGAG")
    tails = encode.(["ATAACTGTGCGTAATTTGTCGGCATAATATCTAAGTTGAAAACTTCTTTACGCGTTAAGT",
                     "GCTGAATAGGACTAAGACCCGCAGCGTGTGTAGGGCCCTGGGACGCCTCGTTCCCATTTA"])
    targets = [vcat(block, tails[1]), vcat(block, tails[2]), vcat(block, uniq)]
    db = SequenceDB(["t1", "t2", "t3"], targets)
    q = vcat(block, uniq)
    block_kmers = Set(NucleotideSearch._unique_kmers(block, 8))
    query_kmers = Set(NucleotideSearch._unique_kmers(q, 8))
    for target in targets[1:2]
        @test intersect(query_kmers, Set(NucleotideSearch._unique_kmers(target, 8))) == block_kmers
    end
    st = NucleotideSearch._SearchState(3)
    count_cands(maxp, mww) = begin
        s = NucleotideSearch._settings(:query, 80, 2, 0.5, 0, 10, 32, mww,
                                       NucleotideSearch._default_scoring(:query), false;
                                       max_postings = maxp)
        c = Tuple{Int32,Int,Int}[]
        NucleotideSearch._collect_candidates!(c, st, NucleotideSearch._unique_kmers(q, 8), db, s, 1)
        sort!(c)
    end
    all3 = count_cands(typemax(Int), 12)
    @test [x[2] for x in all3] == [3, 1, 2]          # all share the block
    capped = count_cands(1, 12)                       # k-mers in > 1 target skipped
    @test [x[2] for x in capped] == [3]
    # Only the block k-mers are skipped (k-mers across the junction are rare too).
    @test -capped[1][1] == length(NucleotideSearch._unique_kmers(q, 8)) -
                           length(NucleotideSearch._unique_kmers(block, 8))
    # Too few rare k-mers: all are used, as without cap.
    @test count_cands(1, 200) == count_cands(typemax(Int), 200)
    @test -count_cands(1, 200)[1][1] == length(NucleotideSearch._unique_kmers(q, 8))
    @test NucleotideSearch._max_postings(0.2, 10_000) == 50_000
    @test NucleotideSearch._max_postings(0.2, 1_000_000) == 200_000
    @test NucleotideSearch._max_postings(1.0, 10) == typemax(Int)
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

@testset "search: better copy away from the seed window" begin
    # The copy with mismatches has the most seeds on one diagonal; the copy
    # with two single-base indels scores more.
    rng = MersenneTwister(41)
    q = String(rand(rng, collect("ACGT"), 300))
    a = collect(q)
    for i in 15:15:300
        a[i] = first(filter(!=(a[i]), collect("ACGT")))
    end
    best = q[1:100] * "A" * q[101:200] * "C" * q[201:end]
    t = String(rand(rng, collect("ACGT"), 150)) * String(a) *
        String(rand(rng, collect("ACGT"), 500)) * best * String(rand(rng, collect("ACGT"), 150))
    full = align(q, t; mode = :query)
    @test full.score == 556
    for kw in ((;), (minwordmatches = 0, maxrejects = 0))
        h = only(search(q, SequenceDB(["r"], [t]); strand = :plus, min_identity = 80,
                        query_cov = 0.9, cigar = true, kw...))
        @test h.score == full.score
        @test (h.target_start, h.target_end) == (full.target_start, full.target_end) == (951, 1252)
        @test h.cigar == full.cigar
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

@testset "search_self: overhangs and asymmetric end gaps" begin
    # iddef 2 ignores end gaps: a shifted overlap of 80 identical bases is a hit
    # whatever the identity threshold
    rng = MersenneTwister(91)
    q = String(rand(rng, collect("ACGT"), 100))
    t = q[21:end] * String(rand(rng, collect("ACGT"), 20))
    for kw in ((;), (minwordmatches = 0, maxrejects = 0))
        res = search_self(SequenceDB(["q", "t"], [q, t]); min_identity = 99, iddef = 2,
                          query_cov = 0.7, kw...)
        @test only(res[1]).identity == only(res[2]).identity == 100.0
        @test only(res[1]).score == align(q, t).score
        @test (only(res[1]).query_start, only(res[1]).target_end) == (21, 80)
    end
    # Different end gap penalties for query and target: each direction is
    # aligned on its own
    q = "ACGTGACCTAGTACCGATGC"
    t = q * "TTTTTTTTTT"
    sc = Scoring(query_end_gap_open = 0, query_end_gap_extend = 0,
                 target_end_gap_open = 20, target_end_gap_extend = 3)
    res = search_self(SequenceDB(["q", "t"], [q, t]); scoring = sc, iddef = 2,
                      min_identity = 0, query_cov = 0, minwordmatches = 0, maxrejects = 0)
    @test only(res[1]).score == align(q, t; scoring = sc).score == 40
    @test only(res[2]).score == align(t, q; scoring = sc).score == -4
end

@testset "search reports and accessors" begin
    rng = MersenneTwister(17)
    refs = [String(rand(rng, collect("ACGT"), 300)) for _ in 1:40]
    db = SequenceDB(["r$i" for i in 1:40], refs)
    @test sequence_ids(db) == ["r$i" for i in 1:40]
    @test sequence_length(db, 3) == 300
    @test_throws BoundsError sequence_length(db, 41)
    queries = [refs[5][51:250], String(rand(rng, collect("ACGT"), 200))]
    # Default limits: the first query stops at maxaccepts = 1, the second
    # (unrelated) runs out of candidates
    reports = SearchReport[]
    res = search(queries, db; maxaccepts = 1, reports)
    @test length(reports) == 2
    @test reports[1].stop === :maxaccepts && reports[1].accepted == 1 && !reports[1].exhaustive
    @test sequence_ids(db)[only(res[1]).target] == "r5"
    @test reports[2].stop === :exhausted && reports[2].accepted == 0
    @test reports[2].examined == reports[2].candidates
    @test !reports[2].exhaustive          # the k-mer filter dropped most targets
    # Every target and strand examined: exhaustive, and the same best hit
    full = search(queries, db; maxaccepts = 0, maxrejects = 0, minwordmatches = 0, reports)
    @test all(r -> r.exhaustive && r.examined == 80 && r.stop === :exhausted, reports)
    @test only(full[1]).target == only(res[1]).target
    # maxrejects ends the search of an unrelated query
    search(queries[2:2], db; minwordmatches = 0, maxrejects = 5, reports)
    @test only(reports).stop === :maxrejects && only(reports).rejected == 5
    # search_self: the query itself is not a candidate
    self_reports = SearchReport[]
    search_self(db; minwordmatches = 0, maxrejects = 0, reports = self_reports)
    @test all(r -> r.candidates == 39 && r.exhaustive, self_reports)
    # Without `reports` the results are the same
    @test search(queries, db; maxaccepts = 1) == res
end

@testset "search reports: limits on the last candidate" begin
    # Coverage is independent of why the search stopped, for both limits,
    # strand counts and execution modes. Truncating hits does not affect it.
    for threads in (1, 2), strand in (:plus, :both), accepted in (true, false)
        query = accepted ? "ACGTACGT" : "AAAAAAAA"
        target = accepted ? query : "CCCCCCCC"
        limit = strand === :both ? 2 : 1
        limits = accepted ? (; maxaccepts = limit, maxrejects = 0) :
                            (; maxaccepts = 0, maxrejects = limit)
        reports = SearchReport[]
        for n in (1, 2)
            refs = SequenceDB(["r$i" for i in 1:n], fill(target, n); k = 3)
            hits = search([query, query], refs; strand, min_identity = 100,
                          minwordmatches = 0, maxhits = 1, reports, threads, limits...)
            @test all(r -> r.stop === (accepted ? :maxaccepts : :maxrejects), reports)
            @test all(r -> r.examined == limit && r.candidates == n * limit, reports)
            @test all(r -> r.exhaustive == (n == 1), reports)
            @test all(h -> length(h) == Int(accepted), hits)
        end
    end
    for threads in (1, 2), accepted in (true, false)
        seqs = accepted ? ["ACGTACGT", "ACGTACGT"] : ["AAAAAAAA", "CCCCCCCC"]
        db = SequenceDB(["a", "b"], seqs; k = 3)
        reports = SearchReport[]
        limits = accepted ? (; maxaccepts = 1, maxrejects = 0) :
                            (; maxaccepts = 0, maxrejects = 1)
        search_self(db; min_identity = 100, minwordmatches = 0, reports, threads, limits...)
        @test all(r -> r.stop === (accepted ? :maxaccepts : :maxrejects), reports)
        @test all(r -> r.examined == r.candidates == 1 && r.exhaustive, reports)
    end
end

@testset "global bands certify the unrestricted optimum" begin
    ns = NucleotideSearch
    q = "CCGTCGAGGCACGTCGATGCACGTCGATGGACGTCGATGCACGTCCATGCTCGTCGTTGCACGTCGATGCACGTCGATGCACGTCGACGCACGTCGTTGC"
    t = "ACGTCGATGCCCGTCGAGGCACGTCGATGCACGTCGATGGACGTCGATGCACGTCCATGCTCGTCGTTGCACGTCGATGCACGTCGATGCACGTCGACGC"
    @test align(q, t).score == 156
    for threads in (1, 4), threshold in (82, 84)
        reports = SearchReport[]
        hits = search(fill(q, 4), SequenceDB(["t"], [t]); mode = :global,
                      strand = :plus, iddef = 1, min_identity = threshold, reports, threads)
        @test all(isempty, hits)
        @test all(r -> r.exhaustive, reports)
        @test all(isempty, search_self(SequenceDB(["q", "t"], [q, t]);
                                      min_identity = threshold, threads))
    end
    rng = MersenneTwister(4473)
    for trial in 1:500
        a = rand(rng, UInt8(0):UInt8(4), rand(rng, 3:35))
        b = rand(rng, UInt8(0):UInt8(4), rand(rng, 3:35))
        sc = Scoring(match = rand(rng, (-3, 1, 2, 4)), mismatch = rand(rng, (-12, -4, 1)),
                     gap_open = rand(rng, (0, 2, 20)), gap_extend = rand(rng, (0, 1, 3)),
                     query_end_gap_open = rand(rng, (0, 2, 20)),
                     query_end_gap_extend = rand(rng, (0, 1, 3)),
                     target_end_gap_open = rand(rng, (0, 2, 20)),
                     target_end_gap_extend = rand(rng, (0, 1, 3)))
        s = ns._settings(:global, rand(rng, (50, 84, 95, 100)), rand(rng, 0:4),
                         rand(rng, (0.0, 0.7, 0.9, 1.0)), 0, 0, 0, 0, sc, false)
        got = ns._global_stats(_Workspace(), a, b, s)
        full = _align(_Workspace(), a, b, -length(a), length(b), s.scoring)
        @test (got.score, got.matches, got.mismatches, got.columns, got.gap_openings,
               got.astart, got.aend, got.bstart, got.bend) ==
              (full.score, full.matches, full.mismatches, full.columns, full.gap_openings,
               full.astart, full.aend, full.bstart, full.bend)
    end
    # Perfect matches can be certified without allocating the full trace.
    a = encode(repeat("ACGT", 100))
    s = ns._settings(:global, 100, 1, 0.9, 0, 0, 0, 0, Scoring(), false)
    ws = _Workspace()
    @test ns._global_stats(ws, a, a, s).score == 800
    @test length(ws.trace) < length(a)^2
end

@testset "negative match scores use a safe score width" begin
    sc = Scoring(match = -999, mismatch = -1, gap_open = 0, gap_extend = 400,
                 terminal_gap_open = 0, terminal_gap_extend = 400)
    a = encode(repeat("A", 50))
    ks = _kernel_scoring(sc, :global)
    @test NucleotideSearch._score_type(50, 50, ks) === Int32
    @test align(a, a; scoring = sc).score == reference_score(a, a, ks) == -40_000
end

@testset "parallel failures join all workers" begin
    entered = Channel{Nothing}(1); release = Channel{Nothing}(1)
    failed = Channel{Nothing}(1); finished = Channel{Nothing}(1)
    parent = Threads.@spawn try
        NucleotideSearch._run_parallel(2, 0, 2) do _, q
            if q == 1
                take!(entered)
                put!(failed, nothing)
                error("injected worker failure")
            else
                put!(entered, nothing)
                take!(release)
                put!(finished, nothing)
            end
        end
        nothing
    catch e
        e
    end
    try
        @test timedwait(() -> isready(failed), 10) == :ok
        @test timedwait(() -> istaskdone(parent), 0.1) == :timed_out
    finally
        put!(release, nothing)
    end
    @test timedwait(() -> istaskdone(parent), 10) == :ok
    @test fetch(parent) isa TaskFailedException
    @test isready(finished)
end

@testset "reciprocal cache has a fixed entry budget" begin
    ns = NucleotideSearch
    sample = _align(_Workspace(), encode("ACGT"), encode("ACGT"), -4, 4,
                    _kernel_scoring(Scoring(), :global))
    for budget in (0, 1, 7, 20)
        cache = ns._PairCache(4; max_entries = budget)
        for key in 1:100
            @test ns._get_or_align!(() -> sample, cache, key) === sample
            @test sum(length, cache.tables) <= budget
        end
        @test sum(length, cache.tables) == budget
    end
    cache = ns._PairCache(4; max_entries = 7)
    Threads.@threads for key in 1:1_000
        ns._get_or_align!(() -> sample, cache, key)
    end
    @test sum(length, cache.tables) <= 7
    cache = ns._PairCache(1; max_entries = 1)
    calls = Ref(0)
    compute() = (calls[] += 1; sample)
    ns._get_or_align!(compute, cache, 1)
    ns._get_or_align!(compute, cache, 1)
    @test calls[] == 1
    @test isempty(only(cache.tables))
    rng = MersenneTwister(19)
    base = randseq(rng, 80)
    db = SequenceDB(["s$i" for i in 1:15], [mutate(rng, base, 0.05) for _ in 1:15])
    expected = search_self(db; cache_entries = 0, threads = 1, maxhits = 0)
    for budget in (0, 1, 10), threads in (1, 4)
        @test search_self(db; cache_entries = budget, threads, maxhits = 0) == expected
    end
    @test_throws ArgumentError search_self(db; cache_entries = -1)
end

@testset "threads: serial and parallel give the same results" begin
    rng = MersenneTwister(23)
    base = String(rand(rng, collect("ACGT"), 400))
    seqs = [join(rand(rng) < 0.03 ? rand(rng, setdiff("ACGT", c)) : c for c in base) for _ in 1:60]
    db = SequenceDB(["s$i" for i in 1:60], seqs)
    ref = search_self(db; maxhits = 0)
    for t in (1, 2, 7, 100)
        @test search_self(db; maxhits = 0, threads = t) == ref
        @test search(seqs[1:10], db; threads = t) == search(seqs[1:10], db)
    end
    @test search(String[], db; threads = 1) == Vector{Hit}[]
    @test_throws ArgumentError search(seqs[1:2], db; threads = 0)
end

end

@testset "Query windows preserve full-kernel tie rules" begin
    q = "ACATAATCAAAGCAGGCGGGTTTCCCGAACTTAGACCACGTCGCCTATGATTGCT"
    target = "T"^60 * "ACATAATCAAAGCAGCCGGGTTTTCCGAACTTAGACCACGTCGCATATGGTTGCT" *
             "C"^120 * "ACATAATCAAAGCGAGGCGGGTTTTCCCGAGACCTTAAGTACCACGTTCTGCCTATGATTGCT" * "A"^60
    scoring = Scoring(match=2, mismatch=0, gap_open=0, gap_extend=1,
                      terminal_gap_open=20, terminal_gap_extend=1)
    db = SequenceDB(["ref"], [target]; k=3)
    full = NucleotideSearch.align(q, target; mode=:query, scoring)
    @test full.score == 102
    @test full.matches == 55
    expected_identity = 100 * 55 / 63
    for threads in (1, 4), threshold in (0, 50, 85, 87, 88, 90, 92)
        hits = search(q, db; strand=:plus, minwordmatches=0, min_identity=threshold,
                      query_cov=0, maxaccepts=0, maxrejects=0, cigar=true, threads, scoring)
        if threshold > expected_identity
            @test isempty(hits)
        else
            @test length(hits) == 1
            hit = only(hits)
            @test hit.identity ≈ expected_identity
            @test (hit.score, hit.matches, hit.mismatches) == (full.score, full.matches, full.mismatches)
            @test (hit.target_start, hit.target_end, hit.cigar) ==
                  (full.target_start, full.target_end, full.cigar)
        end
    end
    # Even a perfect match outside the seed window can win the kernel's tie.
    repeated = q * "C"^120 * q
    repeated_db = SequenceDB(["ref"], [repeated]; k=3)
    for sc in (Scoring(), scoring)
        expected = NucleotideSearch.align(q, repeated; mode=:query, scoring=sc)
        hit = only(search(q, repeated_db; strand=:plus, min_identity=100, query_cov=1,
                          minwordmatches=0, maxaccepts=0, maxrejects=0, cigar=true, scoring=sc))
        @test (hit.score, hit.target_start, hit.target_end, hit.cigar) ==
              (expected.score, expected.target_start, expected.target_end, expected.cigar)
    end
end
