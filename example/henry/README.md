# Henry's Problem - Saltwater Intrusion Benchmark

## Problem Description

Henry's problem is a classical benchmark for density-driven saltwater intrusion in coastal aquifers. It models the steady-state position of the freshwater-saltwater interface in a confined aquifer with:

- Freshwater inflow from the inland side
- Seawater boundary on the coastal side
- Density-driven circulation due to salinity differences

## Standard Parameters (Henry 1964)

| Parameter | Value | Description |
|-----------|-------|-------------|
| L | 2.0 m | Domain length |
| H | 1.0 m | Domain height |
| K | 1.0e-2 m/s | Hydraulic conductivity |
| n | 0.35 | Porosity |
| qf | 6.6e-5 m/s | Freshwater Darcy flux |
| ρf | 1000 kg/m³ | Freshwater density |
| ρs | 1025 kg/m³ | Seawater density |
| Cs | 35 kg/m³ | Seawater salinity |
| Dm | 1.886e-5 m²/s | Molecular diffusion |

## Model Configuration

The frehg model is configured to run Henry's problem with:

### Domain
- NX = 100, NY = 1, NZ = 50
- dx = dy = dz = 0.02 m
- 2D cross-section in x-z plane

### Boundary Conditions

| Boundary | Flow | Salinity |
|----------|------|----------|
| Left (x-) | Fixed flux = 2.64e-8 m³/s/cell | S = 0 (freshwater) |
| Right (x+) | Fixed head = 1.0 m | S = 35 kg/m³ (seawater) |
| Top (z+) | No flow | Zero gradient |
| Bottom (z-) | No flow | Zero gradient |

### Physics
- Baroclinic flow enabled (density from salinity)
- Density ratio: ρ/ρ₀ = 1 + S × 0.000744
- Viscosity ratio: μ₀/μ = 1/(1 + S × 0.0022)

## Running the Simulation

```bash
# From the frehg build directory
./frehg example/henry example/henry/output

# Or with custom paths
./frehg <input_dir> <output_dir>
```

## Expected Results

At steady state, the simulation should show:

1. **Saltwater wedge**: Dense seawater intrudes along the bottom
2. **Mixing zone**: 10-90% salinity contours form a curved interface
3. **Circulation**: Clockwise flow with freshwater flowing over saltwater
4. **Interface position**: 50% isoline tip at approximately x = 1.2-1.4 m

## Plotting Results

```bash
python plot_henry.py <output_dir> <time_step>

# Example:
python plot_henry.py output 50
```

This will generate:
- `henry_salinity_XXXXXX.png`: Salinity distribution with 50% isoline
- `henry_profiles_XXXXXX.png`: Vertical concentration profiles
- `henry_velocity_XXXXXX.png`: Velocity vectors over salinity field

## Model Capabilities for Henry's Problem

| Feature | Status | Notes |
|---------|--------|-------|
| Groundwater flow | ✅ | Richards equation solver |
| Density-driven flow | ✅ | Baroclinic effects enabled |
| Scalar transport | ✅ | Advection-diffusion |
| Dispersion | ✅ | Full tensor (set to 0 for original Henry) |
| Fixed head BC | ✅ | For seawater boundary |
| Fixed flux BC | ✅ | For freshwater inflow |
| Concentration BC | ✅ | For salinity boundaries |

## Verification

Compare results with:
1. Semi-analytical solution (Simpson & Clement 2003)
2. Other numerical codes (SUTRA, SEAWAT, FEFLOW)
3. Laboratory experiments

### Key Metrics
- Position of 50% isoline at bottom (toe)
- Position of 50% isoline at top
- Total salt mass in domain
- Steady-state time

## References

1. Henry, H.R. (1964). Effects of dispersion on salt encroachment in coastal aquifers. USGS Water-Supply Paper 1613-C.

2. Simpson, M.J., & Clement, T.P. (2003). Theoretical analysis of the worthiness of Henry and Elder problems as benchmarks of density-dependent groundwater flow models. Advances in Water Resources, 26(1), 17-31.

3. Segol, G. (1994). Classic Groundwater Simulations: Proving and Improving Numerical Models. Prentice Hall.


