# Changelog

## 0.1.0 (2026-10-19)

First version. License: MIT.

- `SequenceDB`: encoded sequences with an inverted k-mer index; `save_db`
  and `load_db` (versioned little-endian binary format);
  `NucleotideSearch.read_fasta` (not exported).
- `search`: queries against a database, `:query` (fitting, exact reference
  coordinates) or `:global` mode, both strands, coordinates and optional
  CIGAR strings. `Hit` reports `aligned_columns` and the gap openings
  between the first and the last aligned pair.
- `search_self`: all against all, global alignment; the engine of LULU.jl's
  `build_matchlist` (identical match lists).
- `align`: one optimal alignment without band.
- Banded, vectorised affine-gap kernel taken from LULU.jl, with separate
  end-gap penalties for each sequence and alignment coordinates.
