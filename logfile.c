// =============================================================================
// This is the log file explaining modifications to each version of the code
// =============================================================================

// ----- version 1.0 -----
// 2017-06-07
// This is the base version.

// ----- version 1.1 -----
// 2017-06-09
// 1 - Added the wind model
// 2 - Added the depth-dependent bottom CD
// 3 - Enabled time-dependent open boundary scalar

// ----- version 1.2 -----
// 2017-06-16
// 1 - Enabled the cell edge model
// 2 - Enabled spatially distributed initial condition for scalar

// ----- version 2.0 -----
// 2017-07-13
// 1 - Enabled subgrid bathymetry model without wetting and drying

// ----- version 2.1 -----
// 2017-07-23
// 1 - Enabled subgrid bathymetry model with wetting and drying
// BUGS DETECTED:
// 1 - Instability during subgrid wetting drying
// 2 - Non-conservative scalar when using subgrid model

// ----- version 2.2 -----
// 2017-08-05
// 1 - Enabled predictor-corrector updates of the subgrid area to solve the
//     error caused by nonlinearity
// 2 - Enabled subgrid modification of the continuity equation only, but this
//     method seems unstable
// 3 - It remained a problem how to compute the subgrid face area. Using the
//     min area across a grid cell is adopted so far.

// ----- version 2.3 -----
// 2017-08-16
// 1 - Enabled MPI for scalar transport by fixing the bug in scalar limiter
// 2 - Enabled MPI with subgrid model

// ----- version 2.4 -----
// 2017-08-19
// 1 - Fixed an important bug, which did not allow advection in negative
//     directions. Now the 1st upwind advection scheme is complete.

// ----- version 2.5 -----
// 2017-08-21
// 1 - Allowed 3 types of BC (tidal-inflow with tide on YP/YM faces, tidal-tidal)

// ----- version 2.6 -----
// 2017-08-22
// 1 - Added a trigger that forces the program to quit if max CFL is too high.
// 2 - Allowed selectively saving output variables

// ----- version 2.7 -----
// 2017-08-31
// This version passed the wetting drying test on simple staircase bathymetry.
// 1 - Debugged for the scalar limiter to avoid scalar sink terms
// 2 - Removed the small subgrid variable limit to allow more accurate prediction
//     of the wetting drying front
// 3 - Modified the subgrid variable generation algorithm

// ----- version 2.8 -----
// 2017-09-09
// 1 - Fixed the bug that scalar transport was not compatible with MPI.
//     In version 2.7, the scalar transport across different ranks was blocked.

// ----- version 3.0 -----
// 2017-10-05
// 1 - Added a CFL limiter for nonlinear term following FREHD.
// 2 - Added a limiter that removed outflow from dry cells.
// 3 - Added a limiter that limited the smallest subgrid areas and volumes.
//     Otherwise the matrix could go singular.
// 4 - Added a threshold value when searching for subgrid surface index to
//     avoid endless while loop
// 5 - Changed the way the subgrid areas are computed to better capture the
//     surface connectivity, but this needs more testing.

// ----- version 3.1 -----
// 2017-11-24
// This is only a test version because many new settings still requir more tests
// 1 - Debugged for the depth-dependent drag model. In the previous version CD
//     was overlapped by the const CD when computing the thin layer.
// 2 - Removed the bottom index from the subgrid variables. This variable was
//     used to subtract the subgrid volume that is below the edge elevation.
//     But we found that this method did not do what it was supposed to do.
// 3 - Added an automatic blocking identification function in the preprocessing
//     function.
// 4 - Changed the way the subgrid face areas are computed in the preprocessing
//     function. The searching method in the main model was hence simplified.
//     But this requires more testing.
// 5 - Changed the way the bottom curvature is computed in the preprocessing
//     function. Debugged for the curvature-dependent subgrid drag model, but
//     more testings are needed.
// 6 - Enabled saving the subgrid variables and drag CD.
// 7 - Tried to add the evaporation and rainfall model, but not finished yet.
// 8 - The postprocessing scripts are significantly optimized.

// ----- version 3.2 -----
// 2017-12-10
// This is only a test version because many new settings still require more tests
// 1 - Added a subgrid drag Yh model, which is a method to compute bottom drag
//     CD following Volp_2013. But this method seemed underestimates CD. It
//     does not consider anisotropy of CD, so it requires more testing.
// 2 - Added a combined Cv-Yh model, where the Yh is corrected with Cv. But this
//     method requires more testing.
// 3 - Debugged for the subgrid face area at tidal ghost face, but not sure if
//     this modification is compatible with MPI.

// ----- version 3.3 -----
// 2018-01-18
// 1 - Enabled restart from previous model runs.
// 2 - Removed the Wd limiter on wetting and drying for the subgrid model. This
//     limiter causes incorrect inundation area, but without it the subgrid face
//     areas could be wrong. We need to develop a new limiter in the future.

// ----- version 3.4 -----
// 2018-02-18
// 1 - The Wd limiter from version3.3 was debugged to enable interpolation from
//     zero when the depth is small.
// 2 - The direction of the wind stress was debugged, but it needed further
//     validation.
// 3 - In the block-checking function, a wet region was removed if it was too
//     small.
// 4 - When subgrid model is invoked, the depth (used to compute CD) is Modified
//     as the effective depth (sub->V / (dx*dy)).

// ----- version 3.5 -----
// 2018-03-21
// 1 - The wind term was further debugged.
// 2 - Added functions that removes subgrid methods at wetting/drying front, but
//     the result showed that this method was not very effective.
// 3 - Added functions that enabled subgrid modeling with continuity only. This
//     can be achieved by setting->useSubgrid = 2;.

// ----- version 3.6 -----
// 2018-04-11
// 1 - Added a bottom drag model that models the effects of subgrid scale
//     topography. It can be activated by setting->useSubDrag = 3;.

// ----- version 3.7 -----
// 2018-04-17
// 1 - Enabled monitoring and saving volume loss as an output.
// 2 - Added a lower cutoff of the subgrid area/volume ratio to avoid this
//     ratio becomes too small.

// ----- version 3.8 -----
// 2018-05-06
// 1 - Enabled an option of using TVD scheme for scalar advection.
// 2 - Enabled monitoring of volume loss at every time step.
// 3 - When creating the subgrid variables, the bottom variations are modeled using bottom slope (not curvature) because it allows the Cd to both increase and decrease.

// ----- version 3.9 -----
// 2018-06-17
// 1 - The subgrid drag model (setting->useSubDrag = 3) was further improved to include both sidewall drag and bottom drag models.
// 2 - The function to read and allocate subgrid variables was optimized to reduce the memory taken by each thread when modeling in parallel.
// 3 - The subgrid face areas on the YP face was set to zero except where tidal BC was added. This allows the scalar limiter to be used along the tidal boundary.
// ----- version 4.0 -----
// 20180-09-2
// 1 - Removed the restart function. Instead, user can define spatially-variable initial velocities and surface elevation (just like that for scalar)
// 2 - The wind stress term was re-formulated to use edge depth instead of center depth
// 3 - When using tidal-tidal BC, initial tideP and tideM can be different. In this case, the initial surface will linearly change from tideM to tideP.








