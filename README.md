# NucleotideSearch.jl

Pure Julia search and pairwise alignment of nucleotide sequences (amplicons,
ASVs, OTUs), against each other or against a reference database. No external
programs (BLAST, VSEARCH) and no dependencies outside Base.

- **`search`** fits each query into reference sequences (`:query` mode, the
  default). The whole query is aligned, and gaps at its ends (unaligned
  reference bases) are free, so the hit gives **exact start and end
  coordinates of the query on the reference**, also when the ends of the
  query carry mismatches or indels. BLAST, being local, clips such ends.
  Both strands are searched.
- **`search_self`** searches every sequence against all the others with a
  global alignment, the strategy of `vsearch --usearch_global --self`. It is
  the engine of `build_matchlist` in [LULU.jl](../LULU.jl).
- **`SequenceDB`** holds encoded sequences and an inverted k-mer index;
  `save_db`/`load_db` store it in a binary file.
- **`align`** computes one optimal alignment (no band), with coordinates and
  CIGAR string.

The kernel is a banded affine-gap (Gotoh) dynamic programming filled by
anti-diagonals, vectorised (SIMD) with 16-bit scores (32-bit when needed).
Queries run on all threads; results do not depend on the number of threads.

Julia ≥ 1.10. License: MIT. Used by [LULU.jl](../LULU.jl) and
[ASVCoordinates.jl](../ASVCoordinates.jl).

## Installation

The package is not registered. Install it from a local copy:

```julia
using Pkg
Pkg.develop(path = "path/to/NucleotideSearch.jl")
```

## Quick start

```julia
using NucleotideSearch

db = SequenceDB("references.fasta")          # k = 8 by default
save_db("references.nsdb", db)               # later: db = load_db("references.nsdb")

hits = search(["ACGT…", "TTGA…"], db; cigar = true)   # start Julia with -t auto
h = hits[1][1]                                # best hit of the first query
db.ids[h.target], h.strand, h.identity
h.target_start, h.target_end                  # coordinates on the reference
h.query_start, h.query_end                    # < full query: reference too short
h.cigar                                       # e.g. "41S290M"

ml = search_self(SequenceDB(ids, seqs))       # all against all (global)
r = align("ACGTTGCA", "TTACGTAGCATT"; mode = :query)
```

## Search

For each query, targets are ranked by the number of unique k-mers shared with
the query (both strands together, unless `strand = :plus`); targets sharing
fewer than `minwordmatches` are ignored. Candidates are aligned in that
order; a candidate is accepted if its identity is `>= min_identity` and at
least `query_cov` of the query bases are aligned to target bases. The search
stops after `maxaccepts` accepts or `maxrejects` rejects (0: no limit), and
the first `maxhits` hits by decreasing identity are returned.

| Option | `search` | `search_self` |
|---|---|---|
| `mode` | `:query` (or `:global`) | global |
| `min_identity` (percent) | 80 | 84 |
| `iddef` | 2 | 1 |
| `query_cov` | 0.9 | 0.9 |
| `maxhits` / `maxaccepts` / `maxrejects` | 10 / 10 / 32 | 10 / 0 / 32 |
| `minwordmatches` | 12 | 12 |

Identity definitions (`iddef`, as in VSEARCH): `0` matches / shorter length;
`1`, `4` matches / alignment columns; `2` matches / columns excluding end
gaps; `3` `1 - (mismatches + gap openings) / longer length`. Bases other than
A, C, G, T/U never match.

`Scoring` sets match, mismatch and gap penalties (a gap of length `n` costs
`open + n * extend`) and separate penalties for end gaps of the query and of
the target. Defaults: match 2, mismatch −4, gap 20 + 2n, end gaps 2 + n (the
VSEARCH defaults). In `:query` mode the query end gaps are free and the
query ends are **soft-clipped**: an end of the query may be left unaligned,
whether the reference ends there or continues, at a cost of 20 + n by
default (the target end gap penalties). An end is therefore clipped only
after about 7 consecutive mismatches (or about 14 unrelated bases), and a
divergent end (primer remnant, chimera, non-homologous tail) does not distort
the coordinates of the rest. `:global` mode never clips. Only `iddef` 0 and
2 are allowed in `:query` mode (the others count end gaps, which depend on
the window).

In `:query` mode each candidate is aligned in a window of the reference
around the most frequent k-mer diagonal, with a band wide enough for
`min_identity`. If the alignment reaches a window or band boundary inside the
reference, it is recomputed on the whole reference without band.

Coordinates are 1-based and inclusive and refer to the first and last bases
in an aligned pair; query coordinates refer to the query as given, also on
the minus strand. The CIGAR string uses `M`, `I` (query only), `D` (reference
only) and `S` (unaligned query ends).

Notes:

- Near the very ends of the query, shifting the start by a base or two can
  score better than opening a gap (gap open 20): an indel within the first
  or last few bases may move the coordinate by that much. This is a property
  of any affine-gap alignment, not of the band.
- The index takes 4 bytes per unique k-mer of each reference plus the
  sequences (1 byte per base). Measured on synthetic databases of homologous
  650-bp references (100 variants per species, 1–8% divergence), 4 threads:

  | References | Build | Index in memory / file | Load | 500 queries (420 bp) |
  |---|---|---|---|---|
  | 100,000 | 3.3 s | 0.34 GB | 0.2 s | 1.0 s (`maxaccepts = 10`), 4.9 s (100) |
  | 300,000 | 9.5 s | 1.0 GB | 0.5 s | 2.3 s (`maxaccepts = 10`), 6.1 s (100) |

  Memory grows linearly: about 3.3 GB per million 650-bp references, more
  with longer ones. On SILVA v132 (313,502 near full-length 16S, 461 Mbp)
  the index takes 2.4 GB (file 2.3 GB, built in 28 s, loaded in 2 s; peak
  memory 4.9 GB while building). k-mer sampling (minimizers), which would
  reduce memory, is not implemented.
- In large databases with conserved regions (16S: 550 k-mers occur in more
  than half of SILVA, one in 99.7%) every reference becomes a candidate and
  counting shared k-mers dominates the time. `max_kmer_frequency = 0.2`
  skips k-mers found in more than 20% of the references (and in more than
  50,000): on SILVA this halved the search time with the same hits, while
  on databases of 5,000–10,000 references a cap made hits worse, hence the
  50,000 floor. The default of `search` is no cap (`1`); ASVCoordinates.jl
  uses 0.2. Candidates are also ordered by partial sorting, in blocks.
- Bases other than A, C, G, T/U (`N`, IUPAC codes) never match, also in the
  references: an ambiguity code in a reference counts as a mismatch.

## Validation

- `test/runtests.jl` (`julia -t 4 --project=. -e 'using Pkg; Pkg.test()'`,
  about 3,300 tests, 20 s): the kernel against a full-matrix reference
  dynamic programming with random scores, end-gap penalties, soft clipping
  and both score widths, including the score of the traced path;
  coordinates and CIGAR strings on
  synthetic references with mutated query ends, indels, both strands and
  truncated references; windowed search against unbanded alignment;
  database save/load; `search_self`.
- `build_matchlist` of LULU.jl, now built on `search_self`, gives match lists
  identical bit for bit to its previous implementation (example data of the
  R package `lulu`, 2,425 OTUs, nine settings) at the same speed.
- Comparison with BLAST+ 2.17.0 on 303 rbcL ASVs (331 bp) against 5,417
  Diat.barcode references (225–2,874 bp), best BLAST hit of each ASV:

  | | megablast | blastn |
  |---|---|---|
  | best BLAST hits with clipped query ends | 26 | 20 |
  | NucleotideSearch target coordinates equal to BLAST's, after extending BLAST's by the clipped query bases | 297 / 303 | 295 / 303 |

  The remaining cases are references that end inside the ASV: both programs
  give the same coordinates, with the overhang as a soft clip (`41S290M`,
  `326M5S`), except one ASV whose last bases, 3 bases from the end of the
  reference, are ambiguous (BLAST aligns 7 more bases with a gap;
  NucleotideSearch clips 11). `search` of the 303 ASVs against the 5,417
  references took 0.6 s on 4 threads (database built in 0.3 s; loaded in
  0.1 s). See ASVCoordinates.jl for a second data set (COI).

## Origin

The k-mer prefilter and the anti-diagonal kernel come from `build_matchlist`
of LULU.jl (2026), generalised with separate end-gap penalties per sequence,
coordinates, CIGAR strings, a persistent index and the `:query` mode with
soft clipping (a third, vertical gap state for the clipped end, so that the
kernel stays vectorised: 1.6 times the cost of the plain kernel).
