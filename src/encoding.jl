# Nucleotide encoding: A/C/G/T(U) -> 0..3 (any case), anything else -> 4.
# Code 4 never counts as a match and interrupts k-mers.

"""
    encode(s::AbstractString) -> Vector{UInt8}

Encode a nucleotide sequence: `A`, `C`, `G`, `T`/`U` (any case) become
`0`, `1`, `2`, `3`; every other character (`N`, IUPAC codes, gaps) becomes
`4`, which never counts as a match. One code per character.
"""
function encode(s::AbstractString)
    out = Vector{UInt8}(undef, ncodeunits(s))
    n = 0
    for ch in s
        n += 1
        out[n] = ch == 'A' || ch == 'a' ? 0x00 :
                 ch == 'C' || ch == 'c' ? 0x01 :
                 ch == 'G' || ch == 'g' ? 0x02 :
                 ch in ('T', 't', 'U', 'u') ? 0x03 : 0x04
    end
    resize!(out, n)
end

encode(s::AbstractVector{UInt8}) = s

"""
    reverse_complement(s::AbstractVector{UInt8}) -> Vector{UInt8}

Reverse complement of an encoded sequence (code `4` stays `4`).
"""
reverse_complement(s::AbstractVector{UInt8}) =
    UInt8[x < 0x04 ? 0x03 - x : x for x in Iterators.reverse(s)]

# Sorted unique k-mers (2 bits per base) of an encoded sequence; k-mers
# containing code 4 are skipped.
function _unique_kmers(s::AbstractVector{UInt8}, k::Int)
    out = UInt32[]
    mask = (UInt32(1) << (2k)) - UInt32(1)
    w = UInt32(0); valid = 0
    for x in s
        if x > 0x03
            valid = 0; w = UInt32(0)
            continue
        end
        w = ((w << 2) | x) & mask
        valid += 1
        valid >= k && push!(out, w)
    end
    unique!(sort!(out))
end

# K-mers of a sequence with their (1-based) end positions, keeping only k-mers
# that occur once; sorted by k-mer. Used to estimate the alignment diagonal.
function _kmer_positions!(kp::Vector{UInt32}, pos::Vector{Int32}, s::AbstractVector{UInt8}, k::Int)
    empty!(kp); empty!(pos)
    mask = (UInt32(1) << (2k)) - UInt32(1)
    w = UInt32(0); valid = 0
    @inbounds for (p, x) in enumerate(s)
        if x > 0x03
            valid = 0; w = UInt32(0)
            continue
        end
        w = ((w << 2) | x) & mask
        valid += 1
        valid >= k && (push!(kp, w); push!(pos, Int32(p)))
    end
    perm = sortperm(kp; alg = MergeSort)
    kp2 = kp[perm]; pos2 = pos[perm]
    empty!(kp); empty!(pos)
    i = 1; n = length(kp2)
    @inbounds while i <= n
        j = i
        while j < n && kp2[j+1] == kp2[i]
            j += 1
        end
        i == j && (push!(kp, kp2[i]); push!(pos, pos2[i]))
        i = j + 1
    end
    nothing
end
