#!/usr/bin/env python3
"""
Plot results for Henry's Problem - Density-Driven Saltwater Intrusion

This script plots:
1. Salinity distribution (showing the saltwater wedge)
2. Pressure head distribution
3. Time series analysis

Usage:
    python plot_henry.py <output_dir> [time_step]
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import os
import sys
import re

# Domain parameters (should match config.txt)
NX = 100
NY = 1
NZ = 50
dx = 0.02  # m
dz = 0.02  # m
L = NX * dx  # Domain length (2m)
H = NZ * dz  # Domain height (1m)

# Physical parameters
S_sea = 35.0  # Seawater salinity (kg/m³)

def read_vtk_data(filepath):
    """Read data from legacy VTK file format."""
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
        
        # Parse VTK file
        data = {}
        i = 0
        n_data = 0  # Number of data points (either from POINT_DATA or CELL_DATA)
        
        while i < len(lines):
            line = lines[i].strip()
            
            if line.startswith('DIMENSIONS'):
                dims = list(map(int, line.split()[1:]))
                data['dims'] = dims
            elif line.startswith('ORIGIN'):
                data['origin'] = list(map(float, line.split()[1:]))
            elif line.startswith('SPACING'):
                data['spacing'] = list(map(float, line.split()[1:]))
            elif line.startswith('POINT_DATA'):
                n_data = int(line.split()[1])
                data['data_type'] = 'point'
            elif line.startswith('CELL_DATA'):
                n_data = int(line.split()[1])
                data['data_type'] = 'cell'
            elif line.startswith('SCALARS'):
                parts = line.split()
                var_name = parts[1]
                data_type = parts[2] if len(parts) > 2 else 'float'
                i += 1  # Skip LOOKUP_TABLE line
                if i < len(lines) and lines[i].strip().startswith('LOOKUP_TABLE'):
                    i += 1
                # Read the scalar data
                values = []
                while i < len(lines) and len(values) < n_data:
                    val_line = lines[i].strip()
                    if not val_line or val_line.startswith(('SCALARS', 'VECTORS', 'POINT_DATA', 'CELL_DATA')):
                        break
                    vals = val_line.split()
                    try:
                        if data_type == 'int':
                            values.extend([int(v) for v in vals])
                        else:
                            values.extend([float(v) for v in vals])
                    except ValueError:
                        break
                    i += 1
                data[var_name] = np.array(values)
                continue
            
            i += 1
        
        return data
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        import traceback
        traceback.print_exc()
        return None

def reshape_vtk_data(data, var_name):
    """Reshape VTK data to 3D array (NZ, NY, NX)."""
    if var_name not in data:
        return None
    
    values = data[var_name]
    if len(values) == 0:
        return None
        
    dims = data.get('dims', [NX+1, NY+1, NZ+1])
    data_type = data.get('data_type', 'cell')
    
    # VTK structured grids have (nx+1, ny+1, nz+1) nodes for (nx, ny, nz) cells
    # For CELL_DATA, we have cell_dims = dims - 1
    cell_dims = [d - 1 for d in dims]
    n_cells = cell_dims[0] * cell_dims[1] * cell_dims[2]
    n_nodes = dims[0] * dims[1] * dims[2]
    
    try:
        if data_type == 'cell' or len(values) == n_cells:
            # Cell-centered data: reshape to (nz, ny, nx)
            arr = values.reshape((cell_dims[2], cell_dims[1], cell_dims[0]))
        elif len(values) == n_nodes:
            # Node-centered data - take values at cell centers (approx)
            node_arr = values.reshape((dims[2], dims[1], dims[0]))
            arr = node_arr[:-1, :-1, :-1]
        else:
            print(f"Warning: {var_name} size {len(values)} doesn't match cells ({n_cells}) or nodes ({n_nodes})")
            return None
        return arr
    except Exception as e:
        print(f"Warning: Could not reshape {var_name} (size={len(values)}, dims={dims}): {e}")
        return None

def plot_salinity_distribution(salinity, title="Salinity Distribution"):
    """Plot 2D salinity distribution (x-z cross-section)."""
    fig, ax = plt.subplots(figsize=(12, 5))
    
    # Extract x-z plane (y=0)
    if salinity.ndim == 3:
        S = salinity[:, 0, :]  # (NZ, NX)
    else:
        S = salinity.reshape((NZ, NX))
    
    # Create coordinate arrays
    x = np.linspace(dx/2, L - dx/2, NX)
    z = np.linspace(dz/2, H - dz/2, NZ)
    X, Z = np.meshgrid(x, z)
    
    # Plot filled contours
    S_max = max(np.max(S), 0.1)
    levels = np.linspace(0, S_max, 15)
    cf = ax.contourf(X, Z, S, levels=levels, cmap='RdYlBu_r', extend='both')
    
    # Add contour lines for relative concentrations
    if S_max > 0.01:
        cs = ax.contour(X, Z, S, levels=[0.1*S_max, 0.5*S_max, 0.9*S_max], 
                        colors='black', linewidths=1.5)
        ax.clabel(cs, inline=True, fontsize=10, fmt='%.2f')
    
    # Colorbar
    cbar = plt.colorbar(cf, ax=ax, label='Salinity (kg/m³)')
    
    # Labels and title
    ax.set_xlabel('Distance (m)', fontsize=12)
    ax.set_ylabel('Height (m)', fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.set_aspect('equal')
    ax.set_xlim(0, L)
    ax.set_ylim(0, H)
    
    return fig, ax

def plot_pressure_distribution(pressure, title="Pressure Head Distribution"):
    """Plot 2D pressure head distribution."""
    fig, ax = plt.subplots(figsize=(12, 5))
    
    if pressure.ndim == 3:
        P = pressure[:, 0, :]
    else:
        P = pressure.reshape((NZ, NX))
    
    x = np.linspace(dx/2, L - dx/2, NX)
    z = np.linspace(dz/2, H - dz/2, NZ)
    X, Z = np.meshgrid(x, z)
    
    cf = ax.contourf(X, Z, P, levels=20, cmap='viridis')
    cs = ax.contour(X, Z, P, levels=10, colors='white', linewidths=0.5, alpha=0.5)
    
    plt.colorbar(cf, ax=ax, label='Pressure Head (m)')
    
    ax.set_xlabel('Distance (m)', fontsize=12)
    ax.set_ylabel('Height (m)', fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.set_aspect('equal')
    ax.set_xlim(0, L)
    ax.set_ylim(0, H)
    
    return fig, ax

def plot_time_series(output_dir):
    """Plot time series data."""
    ts_file = os.path.join(output_dir, 'time_series.txt')
    if not os.path.exists(ts_file):
        print("Time series file not found")
        return None
    
    # Read time series (skip header lines)
    data = []
    with open(ts_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split()
            if len(parts) >= 4:
                data.append([float(p) for p in parts[:5]])
    
    if not data:
        print("No data in time series file")
        return None
    
    data = np.array(data)
    
    fig, axes = plt.subplots(2, 1, figsize=(10, 8))
    
    # Plot groundwater volume
    axes[0].plot(data[:, 1], data[:, 2], 'b-', linewidth=1.5)
    axes[0].set_xlabel('Time (s)')
    axes[0].set_ylabel('GW Volume (m³)')
    axes[0].set_title('Groundwater Volume over Time')
    axes[0].grid(True, alpha=0.3)
    
    # Plot scalar mass if available
    if data.shape[1] > 4:
        axes[1].plot(data[:, 1], data[:, 4], 'r-', linewidth=1.5)
        axes[1].set_xlabel('Time (s)')
        axes[1].set_ylabel('Scalar Mass (kg)')
        axes[1].set_title('Salt Mass over Time')
        axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    return fig

def main():
    if len(sys.argv) < 2:
        print("Usage: python plot_henry.py <output_dir> [time_step]")
        print("Example: python plot_henry.py out1 80000")
        sys.exit(1)
    
    output_dir = sys.argv[1]
    time_step = int(sys.argv[2]) if len(sys.argv) > 2 else None
    
    print(f"Reading results from: {output_dir}")
    
    # Find available VTK files
    vtk_files = sorted([f for f in os.listdir(output_dir) if f.endswith('.vtk')])
    print(f"Found {len(vtk_files)} VTK files")
    
    if not vtk_files:
        print("No VTK files found!")
        sys.exit(1)
    
    # If no time step specified, use the last one
    if time_step is None:
        # Extract time step from filename (e.g., gw_080000.vtk -> 80000)
        last_file = vtk_files[-1]
        match = re.search(r'_(\d+)\.vtk$', last_file)
        if match:
            time_step = int(match.group(1))
        else:
            time_step = 0
    
    print(f"Plotting time step: {time_step}")
    
    # Find the VTK file for this time step
    vtk_file = None
    for f in vtk_files:
        if f'{time_step:06d}.vtk' in f:
            vtk_file = os.path.join(output_dir, f)
            break
    
    if vtk_file is None:
        print(f"VTK file for time step {time_step} not found")
        print("Available files:")
        for f in vtk_files[:10]:
            print(f"  {f}")
        sys.exit(1)
    
    print(f"Reading: {vtk_file}")
    
    # Read VTK data
    data = read_vtk_data(vtk_file)
    if data is None:
        sys.exit(1)
    
    print(f"Variables found: {[k for k in data.keys() if isinstance(data[k], np.ndarray)]}")
    
    # Plot pressure head (check both 'pressure' and 'pressure_head' names)
    pressure_var = 'pressure_head' if 'pressure_head' in data else 'pressure'
    if pressure_var in data:
        pressure = reshape_vtk_data(data, pressure_var)
        if pressure is not None and pressure.size > 0:
            print(f"Pressure range: {np.min(pressure):.4f} to {np.max(pressure):.4f} m")
            fig1, ax1 = plot_pressure_distribution(pressure, 
                f"Henry's Problem - Pressure Head (step {time_step})")
            fig1.savefig(os.path.join(output_dir, f'henry_pressure_{time_step:06d}.png'), 
                        dpi=150, bbox_inches='tight')
            print(f"Saved: henry_pressure_{time_step:06d}.png")
        else:
            print("Could not read pressure data")
    
    # Plot salinity (scalar_0)
    if 'scalar_0' in data:
        salinity = reshape_vtk_data(data, 'scalar_0')
        if salinity is not None and salinity.size > 0:
            print(f"Salinity range: {np.min(salinity):.4f} to {np.max(salinity):.4f} kg/m³")
            fig2, ax2 = plot_salinity_distribution(salinity, 
                f"Henry's Problem - Salinity (step {time_step})")
            fig2.savefig(os.path.join(output_dir, f'henry_salinity_{time_step:06d}.png'), 
                        dpi=150, bbox_inches='tight')
            print(f"Saved: henry_salinity_{time_step:06d}.png")
        else:
            print("Could not read salinity data")
    
    # Plot time series
    fig3 = plot_time_series(output_dir)
    if fig3:
        fig3.savefig(os.path.join(output_dir, 'henry_time_series.png'),
                    dpi=150, bbox_inches='tight')
        print("Saved: henry_time_series.png")
    
    plt.show()

if __name__ == "__main__":
    main()
