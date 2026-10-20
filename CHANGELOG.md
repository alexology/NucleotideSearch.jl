# Changelog

## Unreleased

- `NucleotideSearch.read_fasta` (and `SequenceDB(fasta_path)`) now uses
  FASTX.jl, the package's first dependency. On SILVA v132 (313,502
  sequences, median of 5 runs) it is 1.4 times faster with one line per
  sequence (0.52 s instead of 0.71 s) and 4.1 times faster with 60-column
  lines (0.74 s instead of 3.07 s), with identical content. All whitespace is now removed from sequence lines (previously only
  at their ends). Malformed files raise an `ArgumentError` naming the file;
  an empty record directly followed by another header, or an empty last
  record without a final newline, is now rejected (FASTX.jl does not accept
  it).

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
