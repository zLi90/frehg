#!/usr/bin/env python3
"""
Convert legacy Frehg benchmark input files (key=value format) to Frehg2 YAML format.

Usage:
    python3 legacy_to_yaml.py --input legacy/benchmarks/b1-sw/input --output benchmarks/b1-sw/b1-sw.yaml
    python3 legacy_to_yaml.py --input legacy/benchmarks/b2-gw/input --output benchmarks/b2-gw/b2-gw.yaml

This is a ONE-TIME conversion for benchmark testing.
Frehg2 does NOT read legacy key=value input files natively.
"""

import sys
import os
import argparse
import re


def parse_legacy_input(filepath):
    """
    Parse a legacy Frehg input file (key=value format).
    
    Legacy format:
        key = value
        key = [v1, v2, ...]
        # comment lines
        key = "string"
    
    Returns a dict of {key: value}.
    """
    params = {}
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            
            # Skip empty lines and comments
            if not line or line.startswith('#'):
                continue
            
            # Parse key = value
            match = re.match(r'^(\w+)\s*=\s*(.+)$', line)
            if not match:
                continue
            
            key = match.group(1).strip()
            value_str = match.group(2).strip()
            
            # Parse value
            if value_str.startswith('"') and value_str.endswith('"'):
                # String
                params[key] = value_str[1:-1]
            elif ',' in value_str:
                # Legacy arrays are usually comma-separated without brackets.
                values = [v.strip() for v in value_str.strip('[]').split(',')]
                try:
                    params[key] = [int(v) for v in values]
                except ValueError:
                    try:
                        params[key] = [float(v) for v in values]
                    except ValueError:
                        params[key] = values
            elif value_str.startswith('[') and value_str.endswith(']'):
                # Array
                inner = value_str[1:-1]
                values = [v.strip() for v in inner.split(',')]
                # Try to convert to numbers
                try:
                    params[key] = [int(v) for v in values]
                except ValueError:
                    try:
                        params[key] = [float(v) for v in values]
                    except ValueError:
                        params[key] = values
            elif '.' in value_str and not value_str.startswith('0x'):
                # Float
                try:
                    params[key] = float(value_str)
                except ValueError:
                    params[key] = value_str
            else:
                # Try int, then float, then string
                try:
                    params[key] = int(value_str)
                except ValueError:
                    try:
                        params[key] = float(value_str)
                    except ValueError:
                        params[key] = value_str
    
    return params


def write_yaml_config(params, output_path, benchmark_type='sw'):
    """
    Write legacy params to YAML format.
    
    benchmark_type: 'sw' (b1-sw) or 'gw' (b2-gw)
    """
    # Build YAML manually to avoid a Python YAML dependency.
    def val(key, default=None):
        return params.get(key, default)

    def bool_text(key, default=0):
        return "true" if int(val(key, default)) != 0 else "false"

    def list_text(key, default):
        value = val(key, default)
        if isinstance(value, list):
            return "[" + ", ".join(str(v) for v in value) + "]"
        return f"[{value}]"

    def flag_or_null(key):
        value = val(key, 0)
        return "null" if value == 0 else str(value)

    lines = ["# Frehg2 Configuration File"]
    lines.append(f"# Converted from legacy input (one-time conversion)")
    lines.append(f"# Benchmark type: {'Surface Water' if benchmark_type == 'sw' else 'Groundwater'}")
    lines.append("")
    
    # --- Simulation section ---
    lines.append("simulation:")
    lines.append(f'  id: "{val("sim_id", "benchmark")}"')
    lines.append(f'  title: "{val("title", "Benchmark")}"')
    lines.append(f'  author: "{val("author", "Frehg2")}"')
    lines.append('  code_version: "2.0.0"')
    lines.append(f'  legacy_input_dir: "{val("finput", "")}"')
    lines.append("")
    
    # --- Domain section ---
    lines.append("domain:")
    lines.append(f'  nx: {val("NX", val("nx", 1))}')
    lines.append(f'  ny: {val("NY", val("ny", 1))}')
    lines.append('  nz: null  # computed by legacy grid rules during initialization')
    lines.append(f'  dx: {val("dx", 1.0)}')
    lines.append(f'  dy: {val("dy", 1.0)}')
    lines.append(f'  dz: {val("dz", 0.1)}')
    lines.append(f'  dz_incre: {val("dz_incre", 1.0)}')
    lines.append(f'  botZ: {val("botZ", 0.0)}')
    lines.append(f'  use_mpi: {bool_text("use_mpi", 0)}')
    lines.append(f'  mpi_nx: {val("mpi_nx", 1)}')
    lines.append(f'  mpi_ny: {val("mpi_ny", 1)}')
    lines.append(f'  bath_file: {bool_text("bath_file", 0)}')
    lines.append(f'  actv_file: {bool_text("actv_file", 0)}')
    lines.append("")
    
    # --- Time section ---
    lines.append("time:")
    lines.append(f'  dt: {val("dt", 1.0)}')
    lines.append(f'  Tend: {val("Tend", 0.0)}')
    lines.append(f'  NT: {val("NT", 0)}')
    lines.append(f'  dt_out: {val("dt_out", 0.0)}')
    lines.append(f'  Co_max: {val("Co_max", 2.0)}')
    lines.append("")
    
    # --- Surface Water section ---
    lines.append("surface_water:")
    lines.append(f'  enable: {bool_text("sim_shallowwater", 1 if benchmark_type == "sw" else 0)}')
    lines.append('  solver: "semi_implicit"')
    lines.append(f'  difuwave: {bool_text("difuwave", 0)}')
    lines.append(f'  init_eta: {val("init_eta", 0.0)}')
    lines.append(f'  eta_file: {flag_or_null("eta_file")}')
    lines.append(f'  uv_file: {flag_or_null("uv_file")}')
    lines.append(f'  bc_type: {list_text("bctype_SW", [0, 0, 0, 0])}  # legacy bctype_SW: x+, x-, y+, y-')
    lines.append(f'  n_tide: {val("n_tide", 0)}')
    lines.append(f'  tide_file: {list_text("tide_file", [])}')
    lines.append(f'  tide_dat_len: {list_text("tide_dat_len", [])}')
    lines.append(f'  tide_locX: {list_text("tide_locX", [])}')
    lines.append(f'  tide_locY: {list_text("tide_locY", [])}')
    lines.append(f'  init_tide: {list_text("init_tide", [])}')
    lines.append(f'  rain_file: {flag_or_null("rain_file")}')
    lines.append(f'  rain_dat_len: {val("rain_dat_len", 0)}')
    lines.append(f'  q_rain: {val("q_rain", 0.0)}')
    lines.append(f'  evap_file: {flag_or_null("evap_file")}')
    lines.append(f'  evap_model: {val("evap_model", 0)}')
    lines.append(f'  q_evap: {val("q_evap", 0.0)}')
    lines.append(f'  evap_dat_len: {val("evap_dat_len", 0)}')
    lines.append(f'  n_inflow: {val("n_inflow", 0)}')
    lines.append(f'  inflow_file: {list_text("inflow_file", [])}')
    lines.append(f'  inflow_dat_len: {list_text("inflow_dat_len", [])}')
    lines.append(f'  inflow_locX: {list_text("inflow_locX", [])}')
    lines.append(f'  inflow_locY: {list_text("inflow_locY", [])}')
    lines.append(f'  init_inflow: {list_text("init_inflow", [])}')
    lines.append(f'  grav: {val("grav", 9.81)}')
    lines.append(f'  viscx: {val("viscx", 1.0e-6)}')
    lines.append(f'  viscy: {val("viscy", 1.0e-6)}')
    lines.append(f'  min_depth: {val("min_dept", 1.0e-8)}')
    lines.append(f'  manning: {val("manning", 0.019)}')
    lines.append(f'  wtfh: {val("wtfh", 1.0e-8)}')
    lines.append(f'  hD: {val("hD", 0.1)}')
    lines.append(f'  rhoa: {val("rhoa", 1.225)}')
    lines.append(f'  rhow: {val("rhow", 998.0)}')
    lines.append(f'  sim_wind: {bool_text("sim_wind", 0)}')
    lines.append(f'  wind_file: {flag_or_null("wind_file")}')
    lines.append(f'  wind_dat_len: {val("wind_dat_len", 0)}')
    lines.append(f'  init_windspd: {val("init_windspd", 0.0)}')
    lines.append(f'  init_winddir: {val("init_winddir", 0.0)}')
    lines.append(f'  Cw: {val("Cw", 0.0013)}')
    lines.append(f'  CwT: {val("CwT", 5.0)}')
    lines.append(f'  north_angle: {val("north_angle", 0.0)}')
    lines.append(f'  use_subgrid: {bool_text("use_subgrid", 0)}')
    lines.append(f'  r_sub: {val("r_sub", 0)}')
    lines.append(f'  eta_sub_min: {val("eta_sub_min", 0.0)}')
    lines.append(f'  eta_sub_max: {val("eta_sub_max", 0.0)}')
    lines.append(f'  deta_sub: {val("deta_sub", 0.0)}')
    lines.append("")
    
    # --- Groundwater section ---
    lines.append("groundwater:")
    lines.append(f'  enable: {bool_text("sim_groundwater", 1 if benchmark_type == "gw" else 0)}')
    lines.append('  solver: "predictor_corrector"')
    lines.append(f'  iter_solve: {val("iter_solve", 0)}')
    lines.append(f'  use_full3d: {bool_text("use_full3d", 0)}')
    lines.append(f'  dt_adjust: {bool_text("dt_adjust", 0)}')
    lines.append(f'  follow_terrain: {bool_text("follow_terrain", 0)}')
    lines.append(f'  sync_coupling: {bool_text("sync_coupling", 1)}')
    lines.append(f'  use_corrector: {bool_text("use_corrector", 1)}')
    lines.append(f'  post_allocate: {bool_text("post_allocate", 0)}')
    lines.append(f'  use_vg: {bool_text("use_vg", 1)}')
    lines.append(f'  use_mvg: {bool_text("use_mvg", 0)}')
    lines.append(f'  aev: {val("aev", -0.02)}')
    lines.append(f'  dt_min: {val("dt_min", 1.0e-4)}')
    lines.append(f'  dt_max: {val("dt_max", 2.0)}')
    lines.append(f'  Co_max: {val("Co_max", 2.0)}')
    lines.append(f'  Ksx: {val("Ksx", 0.0)}')
    lines.append(f'  Ksy: {val("Ksy", 0.0)}')
    lines.append(f'  Ksz: {val("Ksz", 0.0)}')
    lines.append(f'  Ss: {val("Ss", 0.0)}')
    lines.append(f'  soil_a: {val("soil_a", 1.0)}')
    lines.append(f'  soil_n: {val("soil_n", 2.0)}')
    lines.append(f'  wcs: {val("wcs", 0.4)}')
    lines.append(f'  wcr: {val("wcr", 0.0)}')
    lines.append(f'  init_wc: {val("init_wc", 0.0)}')
    lines.append(f'  init_h: {val("init_h", 0.0)}')
    lines.append(f'  init_wt_rel: {val("init_wt_rel", 0.0)}')
    lines.append(f'  init_wt_abs: {val("init_wt_abs", 0.0)}')
    lines.append(f'  h_file: {flag_or_null("h_file")}')
    lines.append(f'  wc_file: {flag_or_null("wc_file")}')
    lines.append(f'  qtop: {val("qtop", 0.0)}')
    lines.append(f'  qbot: {val("qbot", 0.0)}')
    lines.append(f'  htop: {val("htop", 0.0)}')
    lines.append(f'  hbot: {val("hbot", 0.0)}')
    lines.append(f'  qyp: {val("qyp", 0.0)}')
    lines.append(f'  qym: {val("qym", 0.0)}')
    lines.append(f'  bc_type: {list_text("bctype_GW", [0, 0, 0, 0, 0, 0])}  # legacy bctype_GW: x+, x-, y+, y-, bottom(z+), top(z-)')
    lines.append("")
    
    # --- Solute section ---
    lines.append("solute:")
    lines.append(f'  enable: {bool_text("n_scalar", 0)}')
    lines.append(f'  n_scalar: {val("n_scalar", 0)}')
    lines.append(f'  baroclinic: {bool_text("baroclinic", 0)}')
    lines.append(f'  superbee: {bool_text("superbee", 0)}')
    lines.append(f'  scalar_tide_file: {list_text("scalar_tide_file", [])}')
    lines.append(f'  scalar_tide_datlen: {list_text("scalar_tide_datlen", [])}')
    lines.append(f'  scalar_inflow_file: {list_text("scalar_inflow_file", [])}')
    lines.append(f'  scalar_inflow_datlen: {list_text("scalar_inflow_datlen", [])}')
    lines.append(f'  scalar_surf_file: {list_text("scalar_surf_file", [])}')
    lines.append(f'  scalar_subs_file: {list_text("scalar_subs_file", [])}')
    lines.append(f'  init_s_surf: {list_text("init_s_surf", [])}')
    lines.append(f'  init_s_subs: {list_text("init_s_subs", [])}')
    lines.append(f'  s_tide: {list_text("s_tide", [])}')
    lines.append(f'  s_inflow: {list_text("s_inflow", [])}')
    lines.append(f'  s_yp: {list_text("s_yp", [])}')
    lines.append(f'  s_ym: {list_text("s_ym", [])}')
    lines.append(f'  difux: {val("difux", 0.0)}')
    lines.append(f'  difuy: {val("difuy", 0.0)}')
    lines.append(f'  difuz: {val("difuz", 0.0)}')
    lines.append(f'  disp_lon: {val("disp_lon", 0.0)}')
    lines.append(f'  disp_lat: {val("disp_lat", 0.0)}')
    lines.append("")
    
    # --- Output section ---
    lines.append("output:")
    lines.append('  format: "hdf5"')
    lines.append('  filename: "out/output.h5"')
    lines.append(f'  legacy_reference_dir: "{val("foutput", "out/")}"')
    
    if benchmark_type == 'sw':
        lines.append("  variables:")
        lines.append('    - "water_depth"')
        lines.append('    - "water_surface_elevation"')
        lines.append('    - "velocity_x"')
        lines.append('    - "velocity_y"')
    elif benchmark_type == 'gw':
        lines.append("  variables:")
        lines.append('    - "hydraulic_head"')
        lines.append('    - "water_content"')
        lines.append('    - "darcy_flux_x"')
        lines.append('    - "darcy_flux_y"')
        lines.append('    - "darcy_flux_z"')
    
    lines.append("")
    
    # --- Monitor section ---
    lines.append("monitor:")
    lines.append(f'  n_monitor: {val("n_monitor", 0)}')
    n_monitor = int(val("n_monitor", 0) or 0)
    loc_x = val("monitor_locX", [])
    loc_y = val("monitor_locY", [])
    if not isinstance(loc_x, list):
        loc_x = [loc_x]
    if not isinstance(loc_y, list):
        loc_y = [loc_y]
    if n_monitor == 0:
        lines.append("  points: []")
    else:
        lines.append("  points:")
        for idx in range(n_monitor):
            x = loc_x[idx] if idx < len(loc_x) else 0
            y = loc_y[idx] if idx < len(loc_y) else 0
            lines.append(f'    - name: "monitor{idx + 1}"')
            lines.append('      type: "cell"')
            lines.append(f'      locX: {x}')
            lines.append(f'      locY: {y}')
    
    # Write file
    output = '\n'.join(lines)
    with open(output_path, 'w') as f:
        f.write(output)
    
    return output


def main():
    parser = argparse.ArgumentParser(
        description="Convert legacy Frehg input (key=value) to Frehg2 YAML format"
    )
    parser.add_argument("--input", "-i", required=True,
                        help="Path to legacy input file")
    parser.add_argument("--output", "-o", required=True,
                        help="Path to output YAML file")
    parser.add_argument("--type", "-t", choices=["sw", "gw"], default="sw",
                        help="Benchmark type: sw (surface water) or gw (groundwater)")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input):
        print(f"ERROR: Input file not found: {args.input}")
        return 1
    
    print(f"Reading legacy input: {args.input}")
    params = parse_legacy_input(args.input)
    
    print(f"Parsed {len(params)} parameters:")
    for k, v in sorted(params.items()):
        print(f"  {k} = {v}")
    
    print(f"\nWriting YAML output: {args.output}")
    write_yaml_config(params, args.output, args.type)
    
    print("Done!")
    return 0


if __name__ == "__main__":
    sys.exit(main())
