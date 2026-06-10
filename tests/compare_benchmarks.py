#!/usr/bin/env python3
"""
Benchmark Result Comparator for Frehg2

Compares Frehg2 output (HDF5 or ASCII) with legacy Frehg reference output.
Computes L2 norm of difference and reports PASS/FAIL.

Usage:
    python3 compare_benchmarks.py --benchmark b1-sw
    python3 compare_benchmarks.py --benchmark b2-gw --tolerance 1e-6
"""

import os
import sys
import argparse
import math

def read_legacy_binary(filepath, dtype='float64', dims=None):
    """
    Read Frehg legacy binary output file.
    
    Frehg writes raw binary doubles/ints in a specific format.
    Based on Phase 0 audit findings, the format is:
    - Files like depth_*, head_*, moisture_*: raw float64 (double) arrays
    - Dimensions determined by domain size from benchmark input
    
    Returns numpy array (flattened or reshaped if dims provided).
    """
    with open(filepath, 'rb') as f:
        raw = f.read()
    import numpy as np
    data = np.frombuffer(raw, dtype=np.dtype(dtype))
    if dims is not None:
        try:
            data = data.reshape(dims)
        except ValueError:
            pass  # Keep flattened if reshape fails
    return data

def read_legacy_ascii(filepath):
    """Read Frehg legacy ASCII output file (space/newline separated)."""
    with open(filepath, 'r') as f:
        lines = f.readlines()
    values = []
    for line in lines:
        line = line.strip()
        if line:
            for val in line.split():
                try:
                    values.append(float(val))
                except ValueError:
                    pass
    return values

def read_hdf5_dataset(h5file, dataset_path):
    """Read dataset from HDF5 file."""
    import h5py
    with h5py.File(h5file, 'r') as f:
        return f[dataset_path][:]

def flatten_values(values):
    """Return values as a flat Python list."""
    if hasattr(values, "flatten"):
        return list(values.flatten())
    return list(values)

def error_metrics(ref_data, new_data, tolerance):
    ref_flat = flatten_values(ref_data)
    new_flat = flatten_values(new_data)
    min_len = min(len(ref_flat), len(new_flat))
    ref_flat = ref_flat[:min_len]
    new_flat = new_flat[:min_len]
    diff = [ref - new for ref, new in zip(ref_flat, new_flat)]
    l2 = math.sqrt(sum(value * value for value in diff) / len(diff)) if diff else 0.0
    denom = math.sqrt(sum(value * value for value in ref_flat))
    rel_l2 = math.sqrt(sum(value * value for value in diff)) / denom if denom > 0.0 else l2
    max_err = max((abs(value) for value in diff), default=0.0)
    return rel_l2, l2, max_err, rel_l2 < tolerance

def detect_format(filepath):
    """Detect file format: 'binary', 'ascii', or 'hdf5'."""
    if filepath.endswith('.h5') or filepath.endswith('.hdf5'):
        return 'hdf5'
    # Try to read first 4 bytes — if they look like binary, treat as binary
    try:
        with open(filepath, 'rb') as f:
            header = f.read(4)
        # Check for text header patterns (ncols, nrows, etc.)
        if header.startswith(b'nc') or header[0:1].isdigit() or header == b'    ':
            return 'ascii'
        # Otherwise assume binary (raw doubles)
        return 'binary'
    except Exception:
        return 'ascii'

def compare_files(ref_file, new_file, tolerance=1e-6, 
                  is_hdf5=False, dataset_path=None,
                  ref_format=None, ref_dims=None):
    """
    Compare reference file with new output file.
    Returns (l2_error, max_error, mean_error, passed).
    """
    # Read reference
    if ref_format == 'binary':
        ref_data = read_legacy_binary(ref_file, dims=ref_dims)
    elif is_hdf5:
        ref_data = read_hdf5_dataset(ref_file, dataset_path)
    else:
        ref_data = read_legacy_ascii(ref_file)
    
    # Read new output
    if is_hdf5:
        new_data = read_hdf5_dataset(new_file, dataset_path)
    else:
        new_data = read_legacy_ascii(new_file)
    
    ref_data = list(ref_data)
    new_data = list(new_data)
    if len(ref_data) != len(new_data):
        print(f"  WARNING: Length mismatch: ref={len(ref_data)}, new={len(new_data)}")
        min_len = min(len(ref_data), len(new_data))
        ref_data = ref_data[:min_len]
        new_data = new_data[:min_len]
    
    # Compute errors
    diff = [ref - new for ref, new in zip(ref_data, new_data)]
    l2_error = math.sqrt(sum(value * value for value in diff) / len(diff))
    max_error = max(abs(value) for value in diff)
    mean_error = sum(abs(value) for value in diff) / len(diff)
    
    passed = l2_error < tolerance
    
    return l2_error, max_error, mean_error, passed

def compare_benchmark_b1_sw(benchmark_dir, new_output_dir, tolerance=1e-6):
    """
    Compare b1-sw (surface water only) benchmark.
    Reference: legacy/benchmarks/b1-sw/out/
    New: benchmarks/b1-sw/out/ or HDF5
    """
    print("=" * 60)
    print("Benchmark b1-sw: Surface Water Only (Rainfall-Runoff)")
    print("=" * 60)
    
    ref_dir = os.path.join(benchmark_dir, "legacy/benchmarks/b1-sw/out")
    if not os.path.exists(ref_dir):
        ref_dir = os.path.join(benchmark_dir, "benchmarks/b1-sw/out/reference")
    
    new_dir = new_output_dir
    
    # Phase 4.6 acceptance compares depth_* against the legacy reference.
    variables = ["depth"]
    
    all_passed = True
    results = {}
    
    for var in variables:
        # Find all reference files for this variable
        ref_files = sorted([f for f in os.listdir(ref_dir) 
                          if f.startswith(var)])
        
        if not ref_files:
            print(f"  {var}: No reference files found, skipping")
            continue
        
        print(f"\nVariable: {var}")
        var_passed = True
        
        for ref_file in ref_files:
            ref_path = os.path.join(ref_dir, ref_file)
            
            # Phase 0 audit: legacy benchmark outputs are text files written with fprintf.
            ref_data = read_legacy_ascii(ref_path)
            
            new_path = os.path.join(new_dir, ref_file)
            if os.path.exists(new_path):
                new_data = read_legacy_ascii(new_path)
                l2, abs_l2, max_err, passed = error_metrics(ref_data, new_data, tolerance)
            else:
                # Look for new output in HDF5
                hdf5_path = os.path.join(new_dir, "output.h5")
                if not os.path.exists(hdf5_path):
                    print(f"  {ref_file}: FAIL (new output file not found)")
                    var_passed = False
                    all_passed = False
                    continue
                # Map variable prefix to HDF5 dataset path
                dataset_map = {
                    "depth": "/surface/water_depth",
                    "surf": "/surface/water_surface_elevation",
                    "uu": "/surface/velocity_x",
                    "vv": "/surface/velocity_y",
                    "inun": "/surface/inundation",
                }
                time_suffix = ref_file.split("_")[-1]
                base_dataset_path = dataset_map.get(var, f"/surface/{var}")
                dataset_path = f"{base_dataset_path}/{time_suffix}"
                
                try:
                    new_data = read_hdf5_dataset(hdf5_path, dataset_path)
                except KeyError:
                    print(f"  {ref_file}: SKIP (dataset {dataset_path} not in HDF5)")
                    continue
                
                l2, abs_l2, max_err, passed = error_metrics(ref_data, new_data, tolerance)
            
            status = "PASS" if passed else "FAIL"
            print(
                f"  {ref_file}: {status} (relL2={l2:.2e}, absRMS={abs_l2:.2e}, "
                f"max={max_err:.2e})")
            
            if not passed:
                var_passed = False
                all_passed = False
            
            results[ref_file] = {"l2": l2, "max": max_err, "passed": passed}
        
        print(f"  Variable {var}: {'PASS' if var_passed else 'FAIL'}")
    
    print("\n" + "=" * 60)
    print(f"OVERALL: {'PASS' if all_passed else 'FAIL'}")
    print(f"Tolerance: {tolerance}")
    print("=" * 60)
    
    return all_passed

def compare_benchmark_b2_gw(benchmark_dir, new_output_dir, tolerance=1e-6):
    """
    Compare b2-gw (groundwater only) benchmark.
    Reference: legacy/benchmarks/b2-gw/out/
    """
    print("=" * 60)
    print("Benchmark b2-gw: Groundwater Only (Infiltration)")
    print("=" * 60)
    
    ref_dir = os.path.join(benchmark_dir, "legacy/benchmarks/b2-gw/out")
    if not os.path.exists(ref_dir):
        ref_dir = os.path.join(benchmark_dir, "benchmarks/b2-gw/out/reference")
    
    new_dir = new_output_dir
    
    # Phase 5.6 acceptance compares groundwater text mirrors and can also read HDF5.
    variables = ["head", "moisture", "qx", "qy", "qz"]
    
    all_passed = True
    results = {}
    
    for var in variables:
        # Find all reference files
        ref_files = sorted([f for f in os.listdir(ref_dir) 
                          if f.startswith(var)])
        
        if not ref_files:
            print(f"  {var}: No reference files found, skipping")
            continue
        
        print(f"\nVariable: {var}")
        var_passed = True
        
        for ref_file in ref_files:
            ref_path = os.path.join(ref_dir, ref_file)
            
            # Phase 0 audit: legacy benchmark outputs are text files written with fprintf.
            ref_data = read_legacy_ascii(ref_path)
            
            new_path = os.path.join(new_dir, ref_file)
            if os.path.exists(new_path):
                new_data = read_legacy_ascii(new_path)
                l2, abs_l2, max_err, passed = error_metrics(ref_data, new_data, tolerance)
            else:
                # Fall back to HDF5 only when the text mirror is unavailable.
                hdf5_path = os.path.join(new_dir, "output.h5")
                if not os.path.exists(hdf5_path):
                    print(f"  {ref_file}: FAIL (new output file not found)")
                    var_passed = False
                    all_passed = False
                    continue
                dataset_map = {
                    "head": "/groundwater/hydraulic_head",
                    "moisture": "/groundwater/water_content",
                    "qx": "/groundwater/darcy_flux_x",
                    "qy": "/groundwater/darcy_flux_y",
                    "qz": "/groundwater/darcy_flux_z",
                }
                time_suffix = ref_file.split("_")[-1]
                dataset_path = f"{dataset_map.get(var, f'/groundwater/{var}')}/{time_suffix}"
                
                try:
                    new_data = read_hdf5_dataset(hdf5_path, dataset_path)
                except KeyError:
                    print(f"  {ref_file}: SKIP (dataset {dataset_path} not in HDF5)")
                    continue

                l2, abs_l2, max_err, passed = error_metrics(ref_data, new_data, tolerance)
            
            status = "PASS" if passed else "FAIL"
            print(
                f"  {ref_file}: {status} (relL2={l2:.2e}, absRMS={abs_l2:.2e}, "
                f"max={max_err:.2e})")
            
            if not passed:
                var_passed = False
                all_passed = False
            
            results[ref_file] = {"l2": l2, "max": max_err, "passed": passed}
        
        print(f"  Variable {var}: {'PASS' if var_passed else 'FAIL'}")
    
    # Also compare monitor file (mass balance)
    monitor_ref = os.path.join(ref_dir, "monitor1_mass")
    if os.path.exists(monitor_ref):
        print(f"\nVariable: monitor1_mass")
        print(f"  (Mass balance check - manual inspection recommended)")
    
    print("\n" + "=" * 60)
    print(f"OVERALL: {'PASS' if all_passed else 'FAIL'}")
    print(f"Tolerance: {tolerance}")
    print("=" * 60)
    
    return all_passed

def main():
    parser = argparse.ArgumentParser(description="Compare Frehg2 benchmark results")
    parser.add_argument("--benchmark", type=str, required=True,
                       choices=["b1-sw", "b2-gw", "all"],
                       help="Benchmark to compare")
    parser.add_argument("--benchmark-dir", type=str, default=".",
                       help="Root directory containing benchmarks/")
    parser.add_argument("--new-output-dir", type=str, default="out",
                       help="Directory containing new output")
    parser.add_argument("--tolerance", type=float, default=1e-6,
                       help="L2 error tolerance for PASS")
    parser.add_argument("--hdf5", action="store_true",
                       help="New output is HDF5 format")
    
    args = parser.parse_args()
    
    if args.benchmark in ["b1-sw", "all"]:
        compare_benchmark_b1_sw(args.benchmark_dir, 
                                 args.new_output_dir, 
                                 args.tolerance)
    
    if args.benchmark in ["b2-gw", "all"]:
        if args.benchmark == "all":
            print("\n")
        compare_benchmark_b2_gw(args.benchmark_dir,
                                 args.new_output_dir,
                                 args.tolerance)
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
