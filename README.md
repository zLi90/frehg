# frehg
Frehg stands for Fine Resolution Environmental Hydrodynamic and Groundwater model. It is a coupled 2D depth-integrated hydrodynamic and 3D variably-saturated groundwater model.

Currently the Frehg model is still under development. As of 2020-11, the capabilities of Frehg includes:
- Simulate overland flow with either 2D depth-integrated Navier-Stokes equations, or a diffusive wave approximation of the shallow water equation.
- Simulate 3D variably-saturated groundwater flow with the predictor-corrector method.
- Simulate coupled surface-subsurface flow under prescribed tidal elevation, river stage, inflow, rainfall and evaporation rates.
- Simulate transport of passive scalars in both surface and subsurface domains.

Plans for future development include:
- Enable simulation of the subgrid-scale topography.
- Enable more types of boundary conditions for the subsurface domain (currently the side and bottom boundaries must be impermeable).

The hydrodynamic part of Frehg is built based on:<br />
    Li, Z., Hodges, B.R. Model instability and channel connectivity for 2D coastal marsh simulations. Environ Fluid Mech 19, 1309–1338 (2019) doi:10.1007/s10652-018-9623-7 <br />
    Li, Z., Hodges, B.R. Modeling subgrid-scale topographic effects on shallow marsh hydrodynamics and salinity transport. Adv. Water Resour 129, 1-15 (2019) doi:10.1016/j.advwatres.2019.05.004 <br />

The variably-saturated groundwater part of Frehg is built based on:<br />
    Li, Z., Ozgen, I. and Maina, F.Z., (2020), A mass-conservative predictor-corrector solution to the 1D Richards equation with adaptive time control, Journal of Hydrology, under review <br />
    
The following figure validates the frehg model against benchmark problems used by Sulis(2010). The black, red and blue curves are surface runoff of the saturated excess and the infiltration excess examples. The frehg model produces similar results with the CATHY and ParFlow models, with minor deviations likely due to unspecified initial water content in Sulis(2010).
    
![alt text](https://github.com/zLi90/frehg/blob/master/frehg_validation1.png)
    
   

A user manual of frehg should be available in the near future (hopefully...)
