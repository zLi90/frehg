#!/usr/bin/env python3
"""Convert SERGHEI-style benchmark folders to Frehg2 production-schema cases.

The converter is intentionally conservative: it preserves original raster, polygon,
time-series, and reference files; emits a Frehg2 production YAML with current
compatibility fields; and writes a README describing every inferred mapping and
known gap. It is designed for benchmark onboarding, not for silently changing
physics.
"""

from __future__ import annotations

import argparse
import math
import re
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


def strip_comment(line: str) -> str:
    return line.split("//", 1)[0].split("#", 1)[0].strip()


def parse_key_value(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    if not path.exists():
        return values
    for raw in path.read_text().splitlines():
        line = strip_comment(raw)
        if not line:
            continue
        if ":" in line:
            key, value = line.split(":", 1)
        else:
            parts = line.split(maxsplit=1)
            if len(parts) != 2:
                continue
            key, value = parts
        values[key.strip()] = value.strip()
    return values


def scalar(values: dict[str, str], key: str, default: str) -> str:
    return values.get(key, default).split()[0]


def as_float(values: dict[str, str], key: str, default: float) -> float:
    try:
        return float(scalar(values, key, str(default)))
    except ValueError:
        return default


def as_int(values: dict[str, str], key: str, default: int) -> int:
    return int(round(as_float(values, key, float(default))))


def split_number_list(value: str) -> list[float]:
    return [float(item) for item in re.split(r"[;,\s]+", value.strip()) if item]


@dataclass
class RasterInfo:
    ncols: int = 1
    nrows: int = 1
    cellsize: float = 1.0
    nodata: float = -9999.0
    data_start_line: int = 0


def read_raster_info(path: Path) -> RasterInfo:
    info = RasterInfo()
    lines = path.read_text().splitlines()
    header_keys = {
        "ncols": "ncols",
        "n_cols": "ncols",
        "nrows": "nrows",
        "n_rows": "nrows",
        "xllcorner": "origin",
        "yllcorner": "origin",
        "yllcorber": "origin",
        "xllcenter": "origin",
        "yllcenter": "origin",
        "cellsize": "cellsize",
        "nodata_value": "nodata",
    }
    for index, raw in enumerate(lines):
        parts = raw.strip().split()
        if len(parts) < 2:
            continue
        key = parts[0].lower()
        if key not in header_keys:
            info.data_start_line = index
            break
        field_name = header_keys[key]
        value = float(parts[1])
        if field_name == "origin":
            info.data_start_line = index + 1
            continue
        if field_name in {"ncols", "nrows"}:
            setattr(info, field_name, int(value))
        else:
            setattr(info, field_name, value)
        info.data_start_line = index + 1
    return info


def raster_values(path: Path, info: RasterInfo) -> list[float]:
    lines = path.read_text().splitlines()[info.data_start_line :]
    values: list[float] = []
    for line in lines:
        values.extend(float(item) for item in line.split())
    return values


def copy_if_exists(source: Path, destination: Path) -> bool:
    if not source.exists():
        return False
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    return True


def write_esri_raster(path: Path, info: RasterInfo, values: Iterable[float | int]) -> None:
    rows: list[str] = []
    value_list = list(values)
    expected = info.ncols * info.nrows
    if len(value_list) < expected:
        value_list.extend([0] * (expected - len(value_list)))
    for row in range(info.nrows):
        start = row * info.ncols
        rows.append(" ".join(str(value) for value in value_list[start : start + info.ncols]))
    write_text(
        path,
        f"ncols {info.ncols}\n"
        f"nrows {info.nrows}\n"
        "xllcorner 0.0\n"
        "yllcorner 0.0\n"
        f"cellsize {info.cellsize}\n"
        "NODATA_value -9999\n"
        + "\n".join(rows)
        + "\n",
    )


def write_ascii_volume(
    path: Path,
    info: RasterInfo,
    values: Iterable[float | int],
    default_layers: int) -> int:
    value_list = list(values)
    cells_per_layer = info.ncols * info.nrows
    if cells_per_layer <= 0:
        raise ValueError("volume raster dimensions must be positive")
    if len(value_list) % cells_per_layer != 0:
        raise ValueError("volume raster value count does not match horizontal dimensions")
    nlayers = len(value_list) // cells_per_layer
    if nlayers == 0:
        nlayers = default_layers
        value_list.extend([0] * (cells_per_layer * nlayers))

    rows: list[str] = []
    for layer in range(nlayers):
        for row in range(info.nrows):
            start = layer * cells_per_layer + row * info.ncols
            rows.append(" ".join(str(value) for value in value_list[start : start + info.ncols]))

    write_text(
        path,
        f"ncols {info.ncols}\n"
        f"nrows {info.nrows}\n"
        f"nlayers {nlayers}\n"
        "xllcorner 0.0\n"
        "yllcorner 0.0\n"
        f"cellsize {info.cellsize}\n"
        "NODATA_value -9999\n"
        + "\n".join(rows)
        + "\n",
    )
    return nlayers


def soil_id_values(path: Path) -> list[int]:
    values: list[int] = []
    if not path.exists():
        return values
    for raw in path.read_text().splitlines():
        line = strip_comment(raw)
        if not line:
            continue
        parts = line.split()
        if parts[0].lower() == "n_soil":
            continue
        values.extend(int(float(item)) for item in parts)
    return values


def write_column(path: Path, values: Iterable[float]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f"{value:.12g}\n" for value in values))


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def yaml_bool(value: bool) -> str:
    return "true" if value else "false"


@dataclass
class Conversion:
    case_id: str
    source_dir: Path
    output_dir: Path
    input_dir: Path
    modules_surface: bool = False
    modules_groundwater: bool = False
    notes: list[str] = field(default_factory=list)
    gaps: list[str] = field(default_factory=list)

    def out(self, relative: str) -> Path:
        return self.output_dir / relative


@dataclass
class GroundwaterBoundary:
    boundary_id: str
    boundary_type: str
    polygon_file: str
    value: float


def convert_polygon(source: Path, destination: Path) -> bool:
    if not source.exists():
        return False
    points: list[str] = []
    for raw in source.read_text().splitlines():
        line = strip_comment(raw)
        if not line or line.upper().startswith("NPOINTS"):
            continue
        parts = line.split()
        if len(parts) >= 2:
            points.append(f"{parts[0]} {parts[1]}")
    if points and points[0] != points[-1]:
        points.append(points[0])
    write_text(destination, "\n".join(points) + "\n")
    return True


def convert_rainfall(source: Path, destination: Path) -> bool:
    if not source.exists():
        return False
    time_units = "s"
    rain_units = "m/s"
    pairs: list[tuple[float, float]] = []
    for raw in source.read_text().splitlines():
        line = strip_comment(raw)
        if not line:
            continue
        parts = line.split()
        key = parts[0].lower()
        if key == "timeunits" and len(parts) > 1:
            time_units = parts[1].lower()
            continue
        if key == "rainunits" and len(parts) > 1:
            rain_units = parts[1].lower()
            continue
        if key in {"timeseriesnp", "nx", "ny"}:
            continue
        if len(parts) >= 2:
            time_value = float(parts[0])
            rain_value = float(parts[1])
            if time_units in {"h", "hr", "hour", "hours"}:
                time_value *= 3600.0
            if rain_units in {"mm/h", "mm/hr"}:
                rain_value /= 1000.0 * 3600.0
            elif rain_units in {"mm/s"}:
                rain_value /= 1000.0
            pairs.append((time_value, rain_value))
    write_text(destination, "".join(f"{time:.12g} {rate:.12g}\n" for time, rate in pairs))
    return True


def convert_gwbc(source: Path) -> tuple[list[GroundwaterBoundary], list[str]]:
    if not source.exists():
        return [], []

    values = parse_key_value(source)
    warnings: list[str] = []
    if not values:
        return [], warnings

    direction = as_int(values, "direction", 0)
    bctype = as_int(values, "bctype", 0)
    polygon = values.get("polygon", "").split()[0]
    if direction == 6 and bctype == 3 and polygon:
        boundary_id = values.get("id", "groundwater_top_flux").split()[0]
        value = as_float(values, "bcvals", 0.0)
        polygon_path = f"polygons/{Path(polygon).stem}.poly"
        return [
            GroundwaterBoundary(
                boundary_id=boundary_id,
                boundary_type="fixed_flow_rate",
                polygon_file=polygon_path,
                value=value,
            )
        ], warnings

    warnings.append("Detailed gwbc.input semantics are preserved but not fully translated.")
    return [], warnings


def soil_types_from_vg(path: Path) -> list[dict[str, float]]:
    values = parse_key_value(path)
    alpha = split_number_list(values.get("alpha", "1.0"))
    n_values = split_number_list(values.get("n", "2.0"))
    ks_values = split_number_list(values.get("Ks", "1.0e-5"))
    theta_s = split_number_list(values.get("ThetaS", values.get("Phi", "0.4")))
    theta_r = split_number_list(values.get("ThetaR", "0.0"))
    count = max(len(alpha), len(n_values), len(ks_values), len(theta_s), len(theta_r))
    soils: list[dict[str, float]] = []
    for index in range(count):
        soils.append(
            {
                "alpha": alpha[min(index, len(alpha) - 1)],
                "n": n_values[min(index, len(n_values) - 1)],
                "ks": ks_values[min(index, len(ks_values) - 1)],
                "theta_s": theta_s[min(index, len(theta_s) - 1)],
                "theta_r": theta_r[min(index, len(theta_r) - 1)],
            }
        )
    return soils


def write_yaml(conversion: Conversion, raster: RasterInfo, params: dict[str, str], gw: dict[str, str],
               sw: dict[str, str], soils: list[dict[str, float]], groundwater_ic: str,
               groundwater_boundaries: list[GroundwaterBoundary]) -> None:
    surface = conversion.modules_surface
    groundwater = conversion.modules_groundwater
    sim_length = as_float(params, "simLength", 0.0)
    dt = as_float(gw, "dt_init", 1.0) if groundwater else 1.0
    dt_out = as_float(params, "outFreq", max(1.0, sim_length))
    dz = as_float(gw, "dz_base", 0.1)
    nz = as_int(gw, "ndepth", 1) if groundwater else 1
    bot_z = -as_float(gw, "height", dz * nz)
    manning = 0.03
    min_depth = as_float(sw, "dryDepth", 1.0e-8)
    init_eta = 0.0
    initial_mode = sw.get("initialMode", "dry").lower()
    initial_value = as_float(sw, "initialValue", 0.0)
    if surface and initial_mode in {"h", "h+z"}:
        init_eta = initial_value
    first_soil = soils[0] if soils else {
        "alpha": 1.0, "n": 2.0, "theta_s": 0.4, "theta_r": 0.08, "ks": 1.0e-5
    }
    soil_yaml = "\n".join(
        f"""    - id: {index}
      name: "serghei_soil_{index}"
      vg:
        alpha: {soil['alpha']}
        n: {soil['n']}
        theta_s: {soil['theta_s']}
        theta_r: {soil['theta_r']}
        aev: {as_float(gw, 'aev', -0.02)}
      conductivity:
        Ksx: {soil['ks']}
        Ksy: {soil['ks']}
        Ksz: {soil['ks']}"""
        for index, soil in enumerate(soils or [first_soil])
    )
    groundwater_bc_yaml = "  groundwater: []"
    if groundwater_boundaries:
        groundwater_bc_yaml = "  groundwater:\n" + "\n".join(
            f"""    - id: "{boundary.boundary_id}"
      selector:
        type: "polygon"
        file: "{boundary.polygon_file}"
      type: "{boundary.boundary_type}"
      value: {boundary.value}"""
            for boundary in groundwater_boundaries
        )
    runtime_enable = bool(groundwater_boundaries)
    top_bc_type = 2 if groundwater_boundaries else 1
    yaml_text = f"""# Generated by scripts/serghei_to_frehg2.py from {conversion.source_dir}

simulation:
  id: "{conversion.case_id}"
  mode: "{'coupled' if surface and groundwater else 'surface_water' if surface else 'groundwater'}"
  legacy_input_dir: "input/"

modules:
  surface_water: {yaml_bool(surface)}
  groundwater: {yaml_bool(groundwater)}
  solute: false

domain:
  nx: {raster.ncols}
  ny: {raster.nrows}
  nz: {nz}
  dx: {raster.cellsize}
  dy: {raster.cellsize}
  dz: {dz}
  dz_incre: {as_float(gw, 'dz_multiplier', 1.0)}
  botZ: {bot_z}
  use_mpi: false
  mpi_nx: {as_int(params, 'parNx', 1)}
  mpi_ny: {as_int(params, 'parNy', 1)}
  bathymetry:
    source: "ascii_raster"
    file: "rasters/dem.asc"
  active_mask:
    source: "constant"
    value: 1
  bath_file: true
  actv_file: false

time:
  start: 0.0
  dt: {dt}
  Tend: {sim_length}
  NT: 1
  dt_out: {dt_out}
  Co_max: {as_float(params, 'cfl', 0.5)}

surface_water:
  enable: {yaml_bool(surface)}
  solver: "semi_implicit"
  difuwave: false
  init_eta: {init_eta}
  eta_file: null
  uv_file: null
  bc_type: [0, 0, 0, 0]
  n_tide: 0
  tide_file: [0]
  tide_dat_len: [0]
  tide_locX: [0, 0]
  tide_locY: [0, 0]
  init_tide: [0.0]
  rain_file: 1
  rain_dat_len: 0
  q_rain: 0.0
  evap_file: null
  evap_model: 0
  q_evap: 0.0
  evap_dat_len: 0
  n_inflow: 0
  inflow_file: [0]
  inflow_dat_len: [0]
  inflow_locX: [0, 0]
  inflow_locY: [0, 0]
  init_inflow: [0.0]
  grav: 9.81
  viscx: 1.0e-6
  viscy: 1.0e-6
  min_depth: {min_depth}
  manning: {manning}
  wtfh: {min_depth}
  hD: 0.1
  rhoa: 1.225
  rhow: 998.0
  sim_wind: false
  wind_file: null
  wind_dat_len: 0
  init_windspd: 0.0
  init_winddir: 0.0
  Cw: 0.0013
  CwT: 5.0
  north_angle: 0.0
  use_subgrid: false
  r_sub: 30
  eta_sub_min: 0.0
  eta_sub_max: 0.9
  deta_sub: 0.05

groundwater:
  enable: {yaml_bool(groundwater)}
  solver: "predictor_corrector"
  iter_solve: 0
  use_full3d: {yaml_bool(groundwater)}
  dt_adjust: true
  follow_terrain: false
  sync_coupling: false
  async: {yaml_bool(as_int(gw, 'async', 0) != 0)}
  use_corrector: true
  post_allocate: false
  use_vg: true
  use_mvg: false
  aev: {as_float(gw, 'aev', -0.02)}
  dt_min: {as_float(gw, 'dt_init', 0.01)}
  dt_max: {as_float(gw, 'dt_max', 1.0)}
  Co_max: {as_float(params, 'cfl', 0.5)}
  Ksx: {first_soil['ks']}
  Ksy: {first_soil['ks']}
  Ksz: {first_soil['ks']}
  Ss: 1.0e-5
  soil_a: {first_soil['alpha']}
  soil_n: {first_soil['n']}
  wcs: {first_soil['theta_s']}
  wcr: {first_soil['theta_r']}
  init_wc: {first_soil['theta_r']}
  init_h: 1.0
  init_wt_rel: 0.5
  init_wt_abs: -0.75
  h_file: null
  wc_file: null
  qtop: 0.0
  qbot: 0.0
  htop: 0.0
  hbot: 0.0
  qyp: 0.0
  qym: 0.0
  bc_type: [0, 0, 0, 0, 0, {top_bc_type}]

coupling:
  enable: {yaml_bool(surface and groundwater)}
  mode: "{'async' if as_int(gw, 'async', 0) != 0 else 'sync' if surface and groundwater else 'off'}"
  surface_dt: 1.0
  groundwater_dt: {dt}
  min_depth: {min_depth}

runtime:
  enable: {yaml_bool(runtime_enable)}

initial_conditions:
  surface:
    eta:
      source: "constant"
      value: {init_eta}
  groundwater:
    hydraulic_head:
{groundwater_ic}

boundary_conditions:
  surface: []
{groundwater_bc_yaml}

sources:
  surface:
    - id: "serghei-rainfall"
      selector:
        type: "polygon"
        file: "polygons/domain.poly"
      type: "rainfall"
      rate:
        source: "time_series"
        file: "timeseries/rainfall.txt"
  groundwater: []

soil:
  map:
    source: "ascii_raster_3d"
    file: "rasters/soilID.asc"
  types:
{soil_yaml}

solute:
  enable: false
  n_scalar: 0
  baroclinic: false
  superbee: false
  scalar_tide_file: [0]
  scalar_tide_datlen: [0]
  scalar_inflow_file: [0]
  scalar_inflow_datlen: [0]
  scalar_surf_file: [0]
  scalar_subs_file: [0]
  init_s_surf: [0.0]
  init_s_subs: [0.0]
  s_tide: [0.0]
  s_inflow: [0.0]
  s_yp: [0.0]
  s_ym: [0.0]
  difux: 1.0e-10
  difuy: 1.0e-10
  difuz: 1.0e-10
  disp_lon: 0.1
  disp_lat: 0.01

monitoring:
  points: []
  polygons: []
  mass_balance:
    enable: true

output:
  format: "hdf5"
  filename: "out/output.h5"
  text_mirror: true
  variables:
    - "water_depth"
    - "hydraulic_head"
    - "water_content"

validation:
  benchmark: "{conversion.case_id}"
  reference:
    type: "serghei"
    path: "reference"
  tolerances:
    relative_l2: null
    max_abs: null
  variables: []

legacy_compatibility:
  enable: true
  required_by_current_driver: true
  remove_after: "P18"
  fields:
    - "simulation.legacy_input_dir"
    - "surface_water.bc_type"
    - "groundwater.bc_type"
    - "monitor.n_monitor"

monitor:
  n_monitor: 0
  points: []
"""
    write_text(conversion.out(f"{conversion.case_id}.yaml"), yaml_text)


def convert_case(source_dir: Path, output_dir: Path, case_id: str, force: bool = False) -> Conversion:
    if output_dir.exists() and force:
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    input_dir = source_dir / "input"
    if not input_dir.exists() and (source_dir / "input.light").exists():
        input_dir = source_dir / "input.light"
    conversion = Conversion(case_id=case_id, source_dir=source_dir, output_dir=output_dir, input_dir=input_dir)

    params = parse_key_value(input_dir / "parameters.input")
    gw = parse_key_value(input_dir / "subsurface.input")
    sw = parse_key_value(input_dir / "sw.input")
    conversion.modules_surface = (input_dir / "sw.input").exists()
    conversion.modules_groundwater = (input_dir / "subsurface.input").exists()

    dem_source = input_dir / "dem.input"
    if not dem_source.exists():
        raise FileNotFoundError(f"missing DEM raster: {dem_source}")
    raster = read_raster_info(dem_source)
    dem_values = raster_values(dem_source, raster)
    write_esri_raster(conversion.out("rasters/dem.asc"), raster, dem_values)
    write_column(conversion.out("input/bath"), dem_values)
    conversion.notes.append(
        "Mapped dem.input to normalized rasters/dem.asc and compatibility input/bath."
    )

    nz = as_int(gw, "ndepth", 1) if conversion.modules_groundwater else 1
    groundwater_ic = '      source: "legacy_water_table"\n      init_h: 1.0\n      init_wt_rel: 0.5\n      init_wt_abs: -0.75'
    head_source = input_dir / "head.input"
    if head_source.exists():
        head_info = read_raster_info(head_source)
        head_info.cellsize = raster.cellsize
        nlayers = write_ascii_volume(
            conversion.out("rasters/head.asc"),
            head_info,
            raster_values(head_source, head_info),
            nz)
        conversion.notes.append(f"Mapped head.input to {nlayers}-layer rasters/head.asc.")
        groundwater_ic = '      source: "ascii_raster_3d"\n      file: "rasters/head.asc"'
    elif copy_if_exists(input_dir / "wt.input", conversion.out("rasters/water_table.asc")):
        conversion.notes.append("Copied wt.input to rasters/water_table.asc.")
        groundwater_ic = '      source: "ascii_raster"\n      file: "rasters/water_table.asc"'
    else:
        conversion.gaps.append("No head.input or wt.input found; groundwater IC is placeholder metadata.")

    soil_source = input_dir / "soilID.input"
    if soil_source.exists():
        nlayers = write_ascii_volume(
            conversion.out("rasters/soilID.asc"),
            raster,
            soil_id_values(soil_source),
            nz)
        conversion.notes.append(f"Mapped soilID.input to {nlayers}-layer rasters/soilID.asc.")
    else:
        write_text(
            conversion.out("rasters/soilID.asc"),
            f"ncols {raster.ncols}\nnrows {raster.nrows}\nxllcorner 0.0\nyllcorner 0.0\n"
            f"cellsize {raster.cellsize}\nNODATA_value -9999\n"
            + "\n".join(" ".join("0" for _ in range(raster.ncols)) for _ in range(raster.nrows))
            + "\n",
        )
        conversion.gaps.append("No soilID.input found; wrote uniform soilID raster.")

    polygon_candidates = ["polygon.input", "polygontop.input", "polygonbot.input"]
    for polygon_name in polygon_candidates:
        if convert_polygon(input_dir / polygon_name, conversion.out(f"polygons/{Path(polygon_name).stem}.poly")):
            conversion.notes.append(f"Mapped {polygon_name} to polygons/{Path(polygon_name).stem}.poly.")
    if not any((conversion.out(f"polygons/{Path(name).stem}.poly")).exists() for name in polygon_candidates):
        width = raster.ncols * raster.cellsize
        height = raster.nrows * raster.cellsize
        write_text(
            conversion.out("polygons/domain.poly"),
            f"0 0\n{width} 0\n{width} {height}\n0 {height}\n0 0\n",
        )
        conversion.gaps.append("No polygon input found; wrote full-domain polygon.")
    elif not conversion.out("polygons/domain.poly").exists():
        first_polygon = next(conversion.out(f"polygons/{Path(name).stem}.poly")
                             for name in polygon_candidates
                             if conversion.out(f"polygons/{Path(name).stem}.poly").exists())
        shutil.copy2(first_polygon, conversion.out("polygons/domain.poly"))

    if convert_rainfall(input_dir / "rainfall.input", conversion.out("timeseries/rainfall.txt")):
        shutil.copy2(conversion.out("timeseries/rainfall.txt"), conversion.out("input/rain"))
        conversion.notes.append("Mapped rainfall.input to timeseries/rainfall.txt and compatibility input/rain.")
    else:
        write_text(conversion.out("input/rain"), "0 0\n1 0\n")
        conversion.gaps.append("No rainfall.input found; wrote zero-rain compatibility series.")

    if copy_if_exists(input_dir / "extbc.input", conversion.out("source/extbc.input")):
        conversion.notes.append("Copied extbc.input for later detailed BC mapping.")
        conversion.gaps.append("Detailed extbc.input semantics are preserved but not fully translated.")
    groundwater_boundaries, gwbc_gaps = convert_gwbc(input_dir / "gwbc.input")
    if copy_if_exists(input_dir / "gwbc.input", conversion.out("source/gwbc.input")):
        if groundwater_boundaries:
            conversion.notes.append("Translated gwbc.input fixed flux to polygon groundwater BC metadata.")
        else:
            conversion.notes.append("Copied gwbc.input for later detailed GW BC mapping.")
        conversion.gaps.extend(gwbc_gaps)
    if copy_if_exists(input_dir / "roughness.input", conversion.out("rasters/roughness.asc")):
        conversion.notes.append("Copied roughness.input to rasters/roughness.asc.")

    for reference_name in ["ReferenceData", "data", "discharge", "ponding"]:
        reference_source = source_dir / reference_name
        if reference_source.exists():
            shutil.copytree(reference_source, conversion.out("reference") / reference_name, dirs_exist_ok=True)
            conversion.notes.append(f"Copied {reference_name} to reference/{reference_name}.")
    soils = soil_types_from_vg(input_dir / "vg.input") if (input_dir / "vg.input").exists() else []
    if not soils:
        conversion.gaps.append("No vg.input found; wrote default uniform soil parameters.")

    write_yaml(conversion, raster, params, gw, sw, soils, groundwater_ic, groundwater_boundaries)
    write_readme(conversion)
    return conversion


def write_readme(conversion: Conversion) -> None:
    notes = "\n".join(f"- {note}" for note in conversion.notes) or "- No direct mappings recorded."
    gaps = "\n".join(f"- {gap}" for gap in conversion.gaps) or "- No known conversion gaps."
    runtime_text = (
        "The generated YAML enables runtime inputs for translated polygon groundwater BCs. "
        "Phase P18 should decide which converted cases are executable regression tests and which "
        "remain conversion fixtures until missing runtime features are implemented."
        if any("Translated gwbc.input" in note for note in conversion.notes)
        else "The generated YAML keeps `runtime.enable: false` so current b1/b2-style "
             "compatibility execution is not confused with fully validated coupled SERGHEI physics. "
             "Phase P18 should decide which converted cases are executable regression tests and "
             "which remain conversion fixtures until missing runtime features are implemented."
    )
    text = f"""# {conversion.case_id} Frehg2 Conversion

This folder was generated by `scripts/serghei_to_frehg2.py` from:

`{conversion.source_dir}`

## Generated Files

- `{conversion.case_id}.yaml`: Frehg2 production-schema configuration with compatibility fields.
- `rasters/`: copied SERGHEI raster inputs.
- `polygons/`: converted polygon files.
- `timeseries/`: converted forcing time series when available.
- `input/`: compatibility column files for the current transitional driver.
- `reference/`: copied reference data when available.

## Mapping Notes

{notes}

## Known Gaps

{gaps}

{runtime_text}
"""
    write_text(conversion.out("README.md"), text)


def write_missing_case(output_dir: Path, case_id: str, source_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    write_text(
        output_dir / "README.md",
        f"""# {case_id} Frehg2 Conversion

No source folder was found at `{source_dir}` in this workspace.

P17 created this placeholder report so b3-b6 conversion status is explicit. Add the SERGHEI-style
source folder and rerun:

```sh
python3 scripts/serghei_to_frehg2.py --input {source_dir} --output {output_dir} --case-id {case_id} --force
```
""",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, help="SERGHEI-style benchmark folder")
    parser.add_argument("--output", type=Path, help="Frehg2 benchmark output folder")
    parser.add_argument("--case-id", help="Frehg2 case id, e.g. b3-kirkland")
    parser.add_argument("--force", action="store_true", help="overwrite an existing output folder")
    parser.add_argument(
        "--batch-defaults",
        action="store_true",
        help="convert known legacy/benchmarks b3-b6 folders into benchmarks/",
    )
    args = parser.parse_args()

    if args.batch_defaults:
        root = Path(__file__).resolve().parents[1]
        legacy_root = root / "legacy/benchmarks"
        cases: dict[str, Path] = {}
        for prefix in ["b3", "b4", "b5", "b6"]:
            matches = sorted(path for path in legacy_root.iterdir() if path.is_dir() and path.name.startswith(prefix))
            if matches:
                cases[matches[0].name] = matches[0]
            else:
                cases[prefix] = legacy_root / prefix
        for case_id, source in cases.items():
            output = root / "benchmarks" / case_id
            if source.exists():
                convert_case(source, output, case_id, args.force)
            else:
                write_missing_case(output, case_id, source)
        return 0

    if args.input is None or args.output is None or args.case_id is None:
        parser.error("--input, --output, and --case-id are required unless --batch-defaults is used")
    convert_case(args.input, args.output, args.case_id, args.force)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
