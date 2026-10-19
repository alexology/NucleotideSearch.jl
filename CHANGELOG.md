# Changelog

## 0.1.0 (2026-10-19)

First version. License: MIT.

- `SequenceDB`: encoded sequences with an inverted k-mer index; `save_db`
  and `load_db` (versioned little-endian binary format);
  `NucleotideSearch.read_fasta` (not exported).
- `search`: queries against a database, `:query` (fitting with soft
  clipping, exact reference coordinates) or `:global` mode, both strands,
  coordinates and optional CIGAR strings. `Hit` reports `aligned_columns`
  and the gap openings between the first and the last aligned pair. In
  `:query` mode only `iddef` 0 and 2 are allowed. `max_kmer_frequency`
  skips very frequent k-mers when ranking candidates (large databases);
  candidates are ordered by partial sorting in blocks.
- `search_self`: all against all, global alignment; the engine of LULU.jl's
  `build_matchlist` (identical match lists).
- `align`: one optimal alignment without band.
- Banded, vectorised affine-gap kernel taken from LULU.jl, with separate
  end-gap penalties for each sequence, alignment coordinates and optional
  soft clipping of the first sequence (a third vertical gap state; the
  non-clipping path is compiled separately and unchanged).
