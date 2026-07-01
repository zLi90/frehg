#!/usr/bin/env python3
"""Generate docs/yaml_schema_v2.md from the production YAML example fixtures (P17, Task 17.3.6).

This is the SINGLE source for the schema documentation: it walks the committed v2 example
configs (the benchmark YAMLs) and emits a hierarchical reference of every top-level key and
nested block, the example value(s) seen, the inferred type, and a curated units/description
annotation. The frozen P0 field names are authoritative; this tool documents them, it does not
rename them.

Usage:
    python3 tools/gen_yaml_schema_doc.py [--check]

    (no args)   regenerate docs/yaml_schema_v2.md
    --check     fail (exit 1) if the committed doc is stale w.r.t. the fixtures
"""
import argparse
import os
import sys

try:
    import yaml
except ImportError:  # pragma: no cover - environment guard
    sys.stderr.write("gen_yaml_schema_doc: PyYAML is required (pip install pyyaml)\n")
    sys.exit(2)

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Example fixtures walked to build the reference (order matters for "first seen" examples).
FIXTURES = [
    "benchmarks/b1-sw/b1-sw.yaml",
    "benchmarks/b2-gw/b2-gw.yaml",
    "benchmarks/b0-lake/b0-lake.yaml",
]

OUT_PATH = os.path.join(REPO, "docs", "yaml_schema_v2.md")

REQUIRED_TOP_LEVEL = ["schema_version", "simulation", "domain", "time", "modules", "output"]

# Curated units / descriptions keyed by the dotted path (leaf keys). Missing entries render "—".
UNITS = {
    "schema_version": ("—", 'production schema version; must be "2.0"'),
    "simulation.id": ("—", "run identifier (used for output/monitor file names)"),
    "simulation.title": ("—", "free-text description"),
    "simulation.mode": ("—", "surface_water | groundwater | coupled | solute"),
    "domain.nx": ("cells", "interior cells in x (>0)"),
    "domain.ny": ("cells", "interior cells in y (>0)"),
    "domain.nz": ("cells", "vertical layers (>0)"),
    "domain.dx": ("m", "cell size in x"),
    "domain.dy": ("m", "cell size in y"),
    "domain.dz": ("m", "top vertical layer thickness"),
    "domain.dz_incre": ("—", "vertical layer growth factor"),
    "domain.botz": ("m", "domain bottom elevation (datum)"),
    "domain.follow_terrain": ("bool", "terrain-following vertical grid"),
    "domain.mpi.enabled": ("bool", "use an explicit MPI process grid"),
    "domain.mpi.nx": ("ranks", "process-grid size in x"),
    "domain.mpi.ny": ("ranks", "process-grid size in y"),
    "domain.bathymetry.from_file": ("bool", "load bathymetry from a raster/list file"),
    "time.dt": ("s", "base time step"),
    "time.t_end": ("s", "simulation end time"),
    "time.max_steps": ("steps", "step cap (0 = derive from t_end/dt)"),
    "time.output_interval": ("s", "field-output cadence (0 = off)"),
    "time.dt_checkpoint": ("s", "checkpoint cadence (0/null = off)"),
    "time.max_checkpoints": ("count", "retained checkpoint sets"),
    "modules.surface_water": ("bool", "enable the SWE solver"),
    "modules.groundwater": ("bool", "enable the Richards solver"),
    "modules.solute": ("bool", "enable solute transport"),
    "surface_water.gravity": ("m/s^2", "gravitational acceleration"),
    "surface_water.manning": ("s/m^(1/3)", "Manning roughness"),
    "surface_water.min_depth": ("m", "wet/dry threshold"),
    "surface_water.waterfall_depth": ("m", "waterfall correction threshold"),
    "surface_water.h_diffusion_ref": ("m", "reference depth for horizontal diffusion"),
    "groundwater.solver": ("—", "pca (predictor-corrector) is the production path"),
    "groundwater.full_3d": ("bool", "3-D Darcy assembly vs vertical-only"),
    "groundwater.adaptive_dt": ("bool", "non-CFL adaptive time stepping"),
    "groundwater.use_corrector": ("bool", "flux-based theta corrector (PCA)"),
    "groundwater.use_vg": ("bool", "van Genuchten retention"),
    "groundwater.use_mvg": ("bool", "modified van Genuchten"),
    "groundwater.air_entry_value": ("m", "MVG air-entry pressure head"),
    "groundwater.dt_max": ("s", "adaptive dt upper clamp"),
    "groundwater.dt_min": ("s", "adaptive dt lower clamp"),
    "groundwater.co_max": ("—", "adaptive-dt change criterion"),
    "groundwater.specific_storage": ("1/m", "specific storage Ss"),
    "groundwater.bc_type_gw": ("—", "deprecated 6-int legacy BC [x+,x-,y+,y-,z+bottom,z-top]"),
    "coupling.mode": ("—", "sync | sequential | async"),
    "coupling.surface_dt": ("s", "surface-water coupling window"),
    "coupling.groundwater_dt": ("s", "groundwater coupling sub-step"),
    "soil.map.from_file": ("bool", "per-cell soil-class raster"),
    "output.format": ("—", "hdf5 (default backend)"),
    "output.filename": ("path", "output file (relative to the YAML dir)"),
    "output.io_mode": ("—", "serial_gather | parallel_collective | file_per_rank"),
}


def type_name(v):
    if isinstance(v, bool):
        return "bool"
    if isinstance(v, int):
        return "int"
    if isinstance(v, float):
        return "float"
    if isinstance(v, str):
        return "string"
    if isinstance(v, list):
        return "sequence"
    if isinstance(v, dict):
        return "map"
    if v is None:
        return "null"
    return type(v).__name__


def fmt_example(v):
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, float, str)):
        s = str(v)
        return s if len(s) <= 28 else s[:25] + "..."
    if isinstance(v, list):
        return f"[{len(v)} item(s)]"
    if isinstance(v, dict):
        return "{...}"
    if v is None:
        return "null"
    return str(v)


def merge(dst, src):
    """Merge mapping `src` into `dst`, preserving first-seen example values."""
    for k, v in src.items():
        if isinstance(v, dict):
            node = dst.setdefault(k, {})
            if isinstance(node, dict):
                merge(node, v)
        else:
            dst.setdefault(k, v)


def walk(node, prefix, lines, depth):
    for key in node:
        path = f"{prefix}.{key}" if prefix else key
        val = node[key]
        unit, desc = UNITS.get(path, ("—", ""))
        required = "**" if (depth == 0 and key in REQUIRED_TOP_LEVEL) else ""
        indent = "  " * depth
        if isinstance(val, dict):
            lines.append(f"{indent}- {required}`{key}`{required} (map)"
                         + (f" — {desc}" if desc else ""))
            walk(val, path, lines, depth + 1)
        else:
            t = type_name(val)
            ex = fmt_example(val)
            extra = f" — {desc}" if desc else ""
            lines.append(
                f"{indent}- {required}`{key}`{required} ({t}, unit: {unit}, e.g. `{ex}`){extra}")


def build_doc():
    merged = {}
    used = []
    for rel in FIXTURES:
        p = os.path.join(REPO, rel)
        if not os.path.exists(p):
            continue
        with open(p) as f:
            doc = yaml.safe_load(f)
        if isinstance(doc, dict):
            merge(merged, doc)
            used.append(rel)

    out = []
    out.append("# YAML Schema Reference (v2) — AUTHORITATIVE")
    out.append("")
    out.append("> **Generated** by `tools/gen_yaml_schema_doc.py` from the committed example")
    out.append("> fixtures. Do not edit by hand; rerun the generator. This document is the full,")
    out.append("> authoritative reference (Appendix A of `INTEGRATED_PLAN.md` is a quick summary).")
    out.append("")
    out.append("The production schema is `schema_version: \"2.0\"`. The single production driver")
    out.append("(`Orchestrator::initialize()`) refuses any config whose `schema_version` is not")
    out.append("`\"2.0\"` or that is missing a required top-level section. Field names are the P0-")
    out.append("frozen names; the v2 schema is additive only (no `grid.dims`, module-list, or")
    out.append("`io.dir`). Use `tools/migrate_yaml_v1_to_v2` to upgrade a v1/experimental YAML.")
    out.append("")
    out.append("## Required top-level fields")
    out.append("")
    for k in REQUIRED_TOP_LEVEL:
        unit, desc = UNITS.get(k, ("—", ""))
        out.append(f"- **`{k}`**" + (f" — {desc}" if desc else ""))
    out.append("")
    out.append("Modules must include the three booleans "
               "`modules.{surface_water,groundwater,solute}`.")
    out.append("")
    out.append("## Full key reference")
    out.append("")
    out.append("Bold keys are required at the top level. `unit`/example/description are curated;")
    out.append("types and example values are read from the fixtures.")
    out.append("")
    lines = []
    walk(merged, "", lines, 0)
    out.extend(lines)
    out.append("")
    out.append("## Source fixtures")
    out.append("")
    out.append("This reference was generated from:")
    for rel in used:
        out.append(f"- `{rel}`")
    out.append("")
    out.append("## Migration (v1/experimental -> v2)")
    out.append("")
    out.append("`tools/migrate_yaml_v1_to_v2 <in.yaml> [out.yaml]` applies the rename map "
               "(`src/io/YamlMigration.cpp`):")
    out.append("")
    out.append("| v1 / experimental key | frozen v2 key |")
    out.append("|---|---|")
    out.append("| `time.max_step` | `time.max_steps` |")
    out.append("| `time.Tend` / top-level `Tend` | `time.t_end` |")
    out.append("| `time.dt_out` | `time.output_interval` |")
    out.append("| `grid.*` | `domain.*` |")
    out.append("| `domain.bot_z` | `domain.botz` |")
    out.append("| `output.directory` / `io.dir` | `output.filename` |")
    out.append("| BC/source `type: discharge\\|depth\\|critical` | `bc_discharge\\|bc_depth\\|bc_critical` |")
    out.append("| source `type: inflow\\|rainfall\\|well` | `inflow_rate\\|rainfall_rate\\|extraction_well` |")
    out.append("| `soil.uniform: true` (+ flat params) | `soil.types[0]` explicit class |")
    out.append("")
    out.append("Archived v1 inputs live at `legacy/benchmarks/*/input.v1.yaml`.")
    out.append("")
    return "\n".join(out) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="exit 1 if the committed doc is stale")
    args = ap.parse_args()
    doc = build_doc()
    if args.check:
        if not os.path.exists(OUT_PATH):
            sys.stderr.write("gen_yaml_schema_doc: doc missing; run the generator\n")
            return 1
        with open(OUT_PATH) as f:
            current = f.read()
        if current != doc:
            sys.stderr.write("gen_yaml_schema_doc: docs/yaml_schema_v2.md is stale; regenerate\n")
            return 1
        print("gen_yaml_schema_doc: up to date")
        return 0
    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    with open(OUT_PATH, "w") as f:
        f.write(doc)
    print(f"gen_yaml_schema_doc: wrote {OUT_PATH}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
