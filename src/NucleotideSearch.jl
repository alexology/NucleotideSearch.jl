"""
    NucleotideSearch

Pure Julia search and pairwise alignment of nucleotide sequences (amplicons,
ASVs, OTUs) against each other or against a reference database.

- [`SequenceDB`](@ref): encoded sequences with an inverted k-mer index; it
  can be saved to and loaded from disk ([`save_db`](@ref), [`load_db`](@ref)).
- [`search`](@ref): queries against a database. In `:query` mode each query
  is fitted entirely into the target (gaps at the ends of the query are free),
  which gives exact start and end coordinates on the reference.
- [`search_self`](@ref): every sequence of a database against all the others
  (global alignment), the strategy of `vsearch --usearch_global --self`.
- [`align`](@ref): one pairwise alignment, with coordinates and CIGAR string.

The k-mer prefilter and the banded, vectorised affine-gap kernel come from
`build_matchlist` of LULU.jl.
"""
module NucleotideSearch

import FASTX

export Scoring, SequenceDB, Hit, PairAlignment, SearchReport
export encode, reverse_complement, align, search, search_self, save_db, load_db
export nsequences, sequence_ids, sequence_length

include("encoding.jl")
include("kernel.jl")
include("align.jl")
include("database.jl")
include("search.jl")

end # module
