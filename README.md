# frehg
Frehg stands for Fine Resolution Environmental Hydrodynamic and Groundwater model. A 3D variably-saturated groundwater solver is added to the FrehdC hydrodynamic model to enable simulation of surface-subsurface interaction.

Currently the Frehg model is still under development. The subsurface part of the model has not been completely coded and tested.

Documentation of the existing FrehdC model:
    Li, Z., Hodges, B.R. Model instability and channel connectivity for 2D coastal marsh simulations. Environ Fluid Mech 19, 1309–1338 (2019) doi:10.1007/s10652-018-9623-7
    Li, Z., Hodges, B.R. Modeling subgrid-scale topographic effects on shallow marsh hydrodynamics and salinity transport. Adv. Water Resour 129, 1-15 (2019) doi:10.1016/j.advwatres.2019.05.004

The subsurface solver is coded by extending an existing 1D Richards solver:
    Lai, W., & Ogden, F. L. (2015). A mass-conservative finite volume predictor-corrector solution of the 1D Richards’ equation. Journal of Hydrology, 523, 119–127. doi:10.1016/j.jhydrol.2015.01.053
