# Frehg2 I/O Strategy (P3.2.1a / P3.2.1b)

This document records the resolved decisions for the production output layer: the
parallel-write vs. compression trade-off, the visualization sidecar, provenance, and
crash-safe checkpoints. The output path is behind the `OutputWriter` interface
(`include/frehg2/io/OutputWriter.hpp`) so a netCDF/ADIOS2 backend can be added later
without touching solver code.

## 1. Data model: owned cells + global-index map

Solvers never hand the writer a halo-padded array. Each rank passes its **owned physical
cells** (no halo) for a field, and the `IoLayout` provides, once, the **global physical
index** of every owned cell:

- surface index `= gi + gj*nx`
- subsurface index `= gi + gj*nx + k*nx*ny` (i fastest, then j, then k)

This single model makes every write mode reassemble to the *same* global field regardless
of how the domain was decomposed (the decomposition need not be contiguous in global
ordering). On-disk datasets are always in global physical order with no halo.

## 2. Parallel-write vs. compression decision (`output.io_mode`)

Collective parallel-HDF5 writes of gzip-*filtered* datasets are fragile and often slow, so
the mode is **chosen from config, never assumed**. The selected mode is recorded in
`/simulation/metadata@io_mode`.

| `output.io_mode`        | Who writes                         | Compression | When to use |
|-------------------------|------------------------------------|-------------|-------------|
| `serial_gather` (default) | rank 0 gathers, writes one file  | gzip level 6 | small/medium runs; best storage; not scalable |
| `parallel_collective`   | all ranks, shared file, MPI-IO     | disabled (chunked, uncompressed) | large runs; scalable; storage cost higher |
| `file_per_rank`         | each rank writes its own file      | gzip level 6 per file | very large scale; reassembled via the stored per-rank global-index map (a VDS index is the future single-file view) |

Rationale: storage efficiency (gzip) and write scalability (collective MPI-IO) pull in
opposite directions. `serial_gather` optimizes storage and is the simplest correct option;
`parallel_collective` optimizes write bandwidth by dropping filters (this build has
`H5_HAVE_PARALLEL_FILTERED_WRITES`, but we still default parallel writes to uncompressed
for robustness); `file_per_rank` removes shared-file contention entirely and stores the
global-index map so a reader (or a future VDS sidecar) reassembles the global field.

## 3. Visualization sidecar (XDMF)

On `close()`, rank 0 writes an XDMF 3.0 `.xmf` temporal collection referencing the HDF5
surface datasets (`<file>:/surface/<field>/<time>`) so ParaView/VisIt open results directly
with no custom scripts. For `serial_gather`/`parallel_collective` the sidecar references the
single shared file; for `file_per_rank` a VDS-backed single-file view is the planned
follow-up (the sidecar currently references rank 0's file as a smoke-level entry).

## 4. Provenance (reproducibility / traceability)

Every output file embeds, under `/simulation/metadata`: `title`, `version`, `date`,
`frehg2_version`, `git_sha`, `git_dirty`, `config_sha256` (SHA-256 of the resolved config),
`build_type`, `compiler`, `kokkos_backend`, `solver_backend`, `mpi_ranks`, and `io_mode`.
Each field dataset carries a `units` attribute (`m`, `m/s`, `m^3/m^3`). `git_sha`/
`git_dirty`/`build_type`/`compiler` come from a CMake-generated `frehg2/io/build_info.hpp`;
`config_sha256` is computed at runtime from the resolved YAML.

## 5. Crash-safe checkpoints

Checkpoints are written to a **separate file per checkpoint set** with a true atomic commit:
write `…ckpt.<time>.h5.tmp` → `H5Fflush` → `fsync` → `rename(2)` to the final path. A crash
mid-write leaves the previous committed checkpoint fully intact (the final path is only
replaced by the atomic rename; a stray `.tmp` never affects it). Only the most recent
`time.max_checkpoints` files are kept. Each checkpoint stores its own `config_sha256` and
`git_sha`; `readCheckpoint` reports whether the stored config hash matches the running
config so a restart can warn/abort per policy. Checkpoint integration into the time loop is
done by the driver in P7; P3 provides the primitives.

## 6. Local-environment notes

- HDF5 here is parallel-enabled (`H5_HAVE_PARALLEL`), so all three modes are exercised; the
  HDF5 **C** API is used (`H5Cpp.h` is not available locally).
- `python3` + `h5py` 3.11 validate that output files open in the standard ecosystem
  (`tests/io/check_h5py.py`, registered as `test_hdf5_h5py` when h5py is importable).
