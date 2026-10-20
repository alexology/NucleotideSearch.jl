"""
    SequenceDB(ids, sequences; k = 8)
    SequenceDB(records::AbstractVector{<:Pair}; k = 8)
    SequenceDB(fasta_path; k = 8)

Encoded nucleotide sequences with an inverted index of their unique k-mers
(`3 <= k <= 12`). Sequences are stored once, concatenated; the index takes
4 bytes per unique k-mer of each sequence plus `8 * 4^k` bytes. Save it with
[`save_db`](@ref) to avoid rebuilding it.
"""
struct SequenceDB
    ids::Vector{String}
    data::Vector{UInt8}        # encoded sequences, concatenated
    starts::Vector{Int}        # sequence s is data[starts[s]:starts[s+1]-1]
    k::Int
    offsets::Vector{Int}       # postings of k-mer w: offsets[w+1]+1:offsets[w+2]
    postings::Vector{Int32}    # sequence indices, increasing for each k-mer
end

nsequences(db::SequenceDB) = length(db.ids)
Base.length(db::SequenceDB) = nsequences(db)
_seq(db::SequenceDB, s::Integer) = view(db.data, db.starts[s]:(db.starts[s+1] - 1))
_seqlen(db::SequenceDB, s::Integer) = db.starts[s+1] - db.starts[s]

function Base.show(io::IO, db::SequenceDB)
    print(io, "SequenceDB(", nsequences(db), " sequences, ", length(db.data), " bases, k = ",
          db.k, ")")
end

function _check_k(k::Integer)
    3 <= k <= 12 || throw(ArgumentError("k must be between 3 and 12"))
    Int(k)
end

function SequenceDB(ids::AbstractVector{<:AbstractString}, sequences::AbstractVector; k::Integer = 8)
    length(ids) == length(sequences) ||
        throw(DimensionMismatch("ids and sequences must have the same length"))
    length(ids) < typemax(Int32) || throw(ArgumentError("too many sequences"))
    k = _check_k(k)
    starts = Vector{Int}(undef, length(sequences) + 1)
    starts[1] = 1
    data = UInt8[]
    for (s, x) in enumerate(sequences)
        append!(data, encode(x))
        starts[s+1] = length(data) + 1
    end
    _build(String.(ids), data, starts, k)
end

SequenceDB(records::AbstractVector{<:Pair}; k::Integer = 8) =
    SequenceDB([string(first(r)) for r in records], [last(r) for r in records]; k)

function SequenceDB(fasta_path::AbstractString; k::Integer = 8)
    ids, seqs = read_fasta(fasta_path)
    SequenceDB(ids, seqs; k)
end

# Inverted index in CSR form, in two passes over the sequences (counts, then
# postings) so that the k-mers of all sequences are never held at once.
function _build(ids::Vector{String}, data::Vector{UInt8}, starts::Vector{Int}, k::Int)
    n = length(ids)
    nk = 1 << (2k)
    offsets = zeros(Int, nk + 1)
    for s in 1:n, w in _unique_kmers(view(data, starts[s]:(starts[s+1] - 1)), k)
        offsets[w+2] += 1
    end
    cumsum!(offsets, offsets)
    postings = Vector{Int32}(undef, offsets[end])
    fillpos = offsets[1:nk] .+ 1
    for s in 1:n, w in _unique_kmers(view(data, starts[s]:(starts[s+1] - 1)), k)
        postings[fillpos[w+1]] = s
        fillpos[w+1] += 1
    end
    SequenceDB(ids, data, starts, k, offsets, postings)
end

"""
    NucleotideSearch.read_fasta(path) -> (ids, sequences)

FASTA reader based on FASTX.jl (not exported): the id is the header up to the
first whitespace; sequence lines are concatenated with all whitespace removed
and case preserved. Malformed input, including an empty record directly
followed by another header, raises an `ArgumentError`.
"""
function read_fasta(path::AbstractString)
    ids = String[]; seqs = String[]
    reader = FASTX.FASTAReader(open(path); copy = false)
    try
        for rec in reader
            id = FASTX.identifier(rec)
            if isempty(id)   # header such as "> id description"
                words = split(FASTX.description(rec))
                id = isempty(words) ? "" : first(words)
            end
            push!(ids, String(id))
            push!(seqs, _strip_whitespace(FASTX.sequence(String, rec)))
        end
    catch e
        e isa ErrorException || e isa EOFError || rethrow()
        throw(ArgumentError("$path: invalid FASTA ($(sprint(showerror, e)))"))
    finally
        close(reader)
    end
    ids, seqs
end

_isspace_byte(b::UInt8) = b == UInt8(' ') || b == UInt8('\t') || b == UInt8('\r') ||
                          b == UInt8('\v') || b == UInt8('\f')

_strip_whitespace(s::String) =
    any(_isspace_byte, codeunits(s)) ? String(filter(!_isspace_byte, codeunits(s))) : s

const _MAGIC = b"NSDB"
const _FORMAT_VERSION = UInt32(1)

"""
    save_db(path, db::SequenceDB)

Write the database (sequences and k-mer index) to a binary file, read back by
[`load_db`](@ref). The format is little-endian and independent of the Julia
version.
"""
function save_db(path::AbstractString, db::SequenceDB)
    open(path, "w") do io
        write(io, _MAGIC, htol(_FORMAT_VERSION), htol(Int64(db.k)), htol(Int64(nsequences(db))))
        for id in db.ids
            write(io, htol(Int64(ncodeunits(id))), id)
        end
        _write_vec(io, db.data); _write_vec(io, db.starts)
        _write_vec(io, db.offsets); _write_vec(io, db.postings)
    end
    path
end

function _write_vec(io::IO, v::Vector{T}) where {T}
    write(io, htol(Int64(length(v))))
    ENDIAN_BOM == 0x04030201 ? write(io, v) : write(io, htol.(v))
end

function _read_vec(io::IO, ::Type{T}) where {T}
    n = ltoh(read(io, Int64))
    0 <= n || throw(ArgumentError("corrupt database file"))
    v = read!(io, Vector{T}(undef, n))
    ENDIAN_BOM == 0x04030201 ? v : ltoh.(v)
end

"""
    load_db(path) -> SequenceDB

Read a database written by [`save_db`](@ref).
"""
function load_db(path::AbstractString)
    open(path) do io
        read(io, 4) == _MAGIC || throw(ArgumentError("$path is not a NucleotideSearch database"))
        v = ltoh(read(io, UInt32))
        v == _FORMAT_VERSION ||
            throw(ArgumentError("$path: unsupported database format version $v"))
        k = Int(ltoh(read(io, Int64))); n = Int(ltoh(read(io, Int64)))
        ids = Vector{String}(undef, n)
        for s in 1:n
            ids[s] = String(read(io, Int(ltoh(read(io, Int64)))))
        end
        data = _read_vec(io, UInt8); starts = _read_vec(io, Int)
        offsets = _read_vec(io, Int); postings = _read_vec(io, Int32)
        (length(starts) == n + 1 && length(offsets) == (1 << (2k)) + 1 &&
         offsets[end] == length(postings) && starts[end] == length(data) + 1) ||
            throw(ArgumentError("corrupt database file"))
        SequenceDB(ids, data, starts, k, offsets, postings)
    end
end
