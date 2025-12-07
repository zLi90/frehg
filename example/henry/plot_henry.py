#!/usr/bin/env python3
"""
Plot results for Henry's Problem - Density-Driven Saltwater Intrusion

This script plots:
1. Salinity distribution (showing the saltwater wedge)
2. Velocity vectors
3. Comparison with analytical/semi-analytical solution

Usage:
    python plot_henry.py <output_dir> [time_step]
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
import os
import sys

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

def read_scalar_data(filepath):
    """Read 3D scalar data from output file."""
    try:
        data = np.loadtxt(filepath)
        # Reshape to 3D (assuming row-major order: k, j, i)
        if data.size == NX * NY * NZ:
            return data.reshape((NZ, NY, NX))
        else:
            print(f"Warning: Data size {data.size} doesn't match grid {NX*NY*NZ}")
            return data
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return None

def read_velocity_data(output_dir, time_step):
    """Read velocity components."""
    vx = read_scalar_data(os.path.join(output_dir, f"flux_x_{time_step:06d}.txt"))
    vz = read_scalar_data(os.path.join(output_dir, f"flux_z_{time_step:06d}.txt"))
    return vx, vz

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
    levels = np.linspace(0, S_sea, 15)
    cf = ax.contourf(X, Z, S, levels=levels, cmap='RdYlBu_r', extend='both')
    
    # Add contour lines
    cs = ax.contour(X, Z, S, levels=[0.1*S_sea, 0.5*S_sea, 0.9*S_sea], 
                    colors='black', linewidths=1.5)
    ax.clabel(cs, inline=True, fontsize=10, fmt='%.1f')
    
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

def plot_saltwater_interface(salinity, ax=None):
    """Plot the 50% isoline (saltwater-freshwater interface)."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 5))
    
    # Extract x-z plane
    if salinity.ndim == 3:
        S = salinity[:, 0, :]
    else:
        S = salinity.reshape((NZ, NX))
    
    x = np.linspace(dx/2, L - dx/2, NX)
    z = np.linspace(dz/2, H - dz/2, NZ)
    X, Z = np.meshgrid(x, z)
    
    # Plot 50% isoline (mixing zone)
    cs = ax.contour(X, Z, S/S_sea, levels=[0.5], colors='red', linewidths=2)
    ax.clabel(cs, inline=True, fontsize=10, fmt='50%%')
    
    return ax

def plot_henry_analytical(ax, a=0.263):
    """
    Plot semi-analytical solution for Henry's problem.
    
    The analytical solution gives the position of the 50% isoline.
    For the standard Henry problem:
    - a = qf/(K*ε) where qf is freshwater flux, K is conductivity, ε is porosity
    - For standard values: a ≈ 0.263
    """
    # Simplified analytical approximation for 50% isoline
    # This is a rough approximation - exact solution requires numerical integration
    x_tip = L * (1 - 0.5 * np.exp(-2*a))  # Toe position
    
    # Simple exponential profile approximation
    x = np.linspace(x_tip, L, 50)
    z = H * (1 - np.exp(-(x - x_tip)/(L - x_tip + 0.01)))
    
    ax.plot(x, z, 'k--', linewidth=2, label='Semi-analytical (approx.)')
    ax.legend()
    
    return ax

def plot_velocity_field(vx, vz, salinity, subsample=5):
    """Plot velocity vectors over salinity field."""
    fig, ax = plt.subplots(figsize=(12, 5))
    
    # Plot salinity background
    if salinity.ndim == 3:
        S = salinity[:, 0, :]
    else:
        S = salinity.reshape((NZ, NX))
    
    x = np.linspace(dx/2, L - dx/2, NX)
    z = np.linspace(dz/2, H - dz/2, NZ)
    X, Z = np.meshgrid(x, z)
    
    cf = ax.contourf(X, Z, S, levels=15, cmap='RdYlBu_r', alpha=0.5)
    plt.colorbar(cf, ax=ax, label='Salinity (kg/m³)')
    
    # Extract velocity components
    if vx is not None and vz is not None:
        if vx.ndim == 3:
            Vx = vx[:, 0, :]
            Vz = vz[:, 0, :]
        else:
            Vx = vx.reshape((NZ, NX))
            Vz = vz.reshape((NZ, NX))
        
        # Subsample for clearer visualization
        X_sub = X[::subsample, ::subsample]
        Z_sub = Z[::subsample, ::subsample]
        Vx_sub = Vx[::subsample, ::subsample]
        Vz_sub = Vz[::subsample, ::subsample]
        
        # Normalize velocity magnitude
        V_mag = np.sqrt(Vx_sub**2 + Vz_sub**2)
        V_max = np.max(V_mag)
        if V_max > 0:
            Vx_norm = Vx_sub / V_max
            Vz_norm = Vz_sub / V_max
        else:
            Vx_norm = Vx_sub
            Vz_norm = Vz_sub
        
        ax.quiver(X_sub, Z_sub, Vx_norm, Vz_norm, color='black', alpha=0.7)
    
    ax.set_xlabel('Distance (m)', fontsize=12)
    ax.set_ylabel('Height (m)', fontsize=12)
    ax.set_title("Velocity Field over Salinity Distribution", fontsize=14)
    ax.set_aspect('equal')
    ax.set_xlim(0, L)
    ax.set_ylim(0, H)
    
    return fig, ax

def plot_concentration_profiles(salinity, x_positions=[0.5, 1.0, 1.5, 1.9]):
    """Plot vertical concentration profiles at different x positions."""
    fig, ax = plt.subplots(figsize=(8, 6))
    
    if salinity.ndim == 3:
        S = salinity[:, 0, :]
    else:
        S = salinity.reshape((NZ, NX))
    
    z = np.linspace(dz/2, H - dz/2, NZ)
    x = np.linspace(dx/2, L - dx/2, NX)
    
    colors = plt.cm.viridis(np.linspace(0, 1, len(x_positions)))
    
    for x_pos, color in zip(x_positions, colors):
        # Find nearest grid index
        i = int(x_pos / dx)
        if i >= NX:
            i = NX - 1
        
        ax.plot(S[:, i]/S_sea, z, color=color, linewidth=2, 
                label=f'x = {x_pos:.1f} m')
    
    ax.set_xlabel('Relative Salinity (S/S_sea)', fontsize=12)
    ax.set_ylabel('Height (m)', fontsize=12)
    ax.set_title('Vertical Salinity Profiles', fontsize=14)
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 1.1)
    ax.set_ylim(0, H)
    
    return fig, ax

def main():
    if len(sys.argv) < 2:
        print("Usage: python plot_henry.py <output_dir> [time_step]")
        print("Example: python plot_henry.py output 50")
        sys.exit(1)
    
    output_dir = sys.argv[1]
    time_step = int(sys.argv[2]) if len(sys.argv) > 2 else 50
    
    print(f"Reading results from: {output_dir}")
    print(f"Time step: {time_step}")
    
    # Read salinity data
    salinity_file = os.path.join(output_dir, f"scalar_0_{time_step:06d}.txt")
    if not os.path.exists(salinity_file):
        # Try alternative naming
        salinity_file = os.path.join(output_dir, f"gw_scalar_0_{time_step:06d}.txt")
    
    if os.path.exists(salinity_file):
        salinity = read_scalar_data(salinity_file)
        if salinity is not None:
            print(f"Salinity range: {np.min(salinity):.2f} to {np.max(salinity):.2f} kg/m³")
            
            # Plot salinity distribution
            fig1, ax1 = plot_salinity_distribution(salinity, 
                f"Henry's Problem - Saltwater Intrusion (t = {time_step * 1000}s)")
            plot_saltwater_interface(salinity, ax1)
            fig1.savefig(os.path.join(output_dir, f'henry_salinity_{time_step:06d}.png'), 
                        dpi=150, bbox_inches='tight')
            print(f"Saved: henry_salinity_{time_step:06d}.png")
            
            # Plot concentration profiles
            fig2, ax2 = plot_concentration_profiles(salinity)
            fig2.savefig(os.path.join(output_dir, f'henry_profiles_{time_step:06d}.png'),
                        dpi=150, bbox_inches='tight')
            print(f"Saved: henry_profiles_{time_step:06d}.png")
            
            # Try to plot velocity field
            vx, vz = read_velocity_data(output_dir, time_step)
            if vx is not None:
                fig3, ax3 = plot_velocity_field(vx, vz, salinity)
                fig3.savefig(os.path.join(output_dir, f'henry_velocity_{time_step:06d}.png'),
                            dpi=150, bbox_inches='tight')
                print(f"Saved: henry_velocity_{time_step:06d}.png")
            
            plt.show()
    else:
        print(f"Error: Salinity file not found: {salinity_file}")
        print("Available files in output directory:")
        if os.path.exists(output_dir):
            for f in sorted(os.listdir(output_dir))[:20]:
                print(f"  {f}")
        sys.exit(1)

if __name__ == "__main__":
    main()

