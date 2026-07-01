# Migrating from legacy `frehg` to Frehg2

The legacy `frehg` model is [deprecated](../legacy/frehg/DEPRECATED.md). This guide shows how to
move a legacy run to Frehg2.

## 1. Input files: legacy `input` → Frehg2 v2 YAML

The legacy code reads a single `key = value` `input` file plus referenced ASCII inputs
(`bath`, `rain`, …). Frehg2 reads a structured **v2 YAML** config (schema in
[`docs/yaml_schema_v2.md`](yaml_schema_v2.md); the legacy→YAML field mapping was frozen in
[`docs/legacy_audit/yaml_schema.md`](legacy_audit/yaml_schema.md)).

Convert automatically with the P0 tool:

```bash
python3 scripts/legacy_to_yaml.py legacy/benchmarks/b1-sw/input benchmarks/b1-sw \
    --id b1-sw --finput-dir legacy/benchmarks/b1-sw/b1-input
# b2-gw additionally takes --warrick-ref <warrick_csv>
```

This writes `<out_dir>/<id>.yaml` in the frozen schema (including a `legacy_raw` block that
preserves *every* original key for 100% round-trip provenance) and copies the referenced ASCII
inputs (`bath`, `rain`, and for `b2-gw` the Warrick reference CSV).

The six shipped benchmark configs in `benchmarks/` were produced this way (b1/b2 in P0; the
SERGHEI/legacy b3–b6 in P18). Use them as worked examples.

### If you already have a v1 YAML

Older Frehg2 v1 YAML configs upgrade to v2 with the migration tool (idempotent):

```bash
./build/tools/migrate_yaml_v1_to_v2 old_config.v1.yaml new_config.yaml
```

The Orchestrator enforces `schema_version: '2.0'` at load (`validateSchemaV2`), so v1 files must be
migrated before they run.

## 2. Running

Legacy:

```bash
make -C legacy/frehg          # archival opt-in only; produces legacy/frehg/frehg
./legacy/frehg/frehg          # reads ./input
```

Frehg2:

```bash
cmake -S . -B build -DCMAKE_C_COMPILER=gcc-15 -DCMAKE_CXX_COMPILER=g++-15 -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/frehg2 benchmarks/b1-sw/b1-sw.yaml
# multi-rank (local MPICH): /Users/zhili/Codes/local/bin/mpiexec -n 2 ./build/frehg2 <config.yaml>
```

## 3. Outputs

| Aspect            | Legacy                                  | Frehg2                                                        |
|-------------------|-----------------------------------------|--------------------------------------------------------------|
| Field output      | per-variable ASCII text snapshots       | HDF5 (`/surface/*`, `/subsurface/*`) + XDMF sidecar (ParaView/VisIt) |
| Provenance        | none                                    | `git_sha`, `config_sha256`, units, build/backend in every file |
| Run summary       | stdout                                  | `simulation_summary.txt` (modules, steps, water mass balance) |
| Monitors          | hard-coded probe prints                 | CSV at `output_interval` (`monitors.*` / `monitoring.points`) |
| Checkpoint/restart| none                                    | atomic HDF5 checkpoints; `--restart <file> --restart-time <t>` |

Legacy text snapshots are deprecated. To compare a Frehg2 run against a legacy reference, use the
validation harness (`tools/run_validation.py`) or `scripts/compare_with_legacy.py`.

## 4. Behavioural notes / parity caveats

- **b1-sw / b2-gw** reproduce the legacy results to the gated tolerances (b1-sw relative depth
  L2 = 0 vs the committed reference; b2-gw element-wise legacy parity via ctest `test_re_b2_gw`).
- **b3–b6** are **review-tier** ports: the current Frehg2 solver scope does not reproduce every
  legacy/SERGHEI physics option (lateral 2-D Richards, Chezy friction, Newton/baroclinic transport).
  See [`docs/benchmarks/`](benchmarks/) for the per-benchmark approximations and the three-tier
  thresholds, and [`docs/validation.md`](validation.md) for the harness.
- Indexing differs: legacy interior ordering vs Frehg2's halo-padded flat layout. Never compare raw
  indices — use `include/frehg2/core/LegacyIndexAdapter.hpp` (the comparison scripts already do).
