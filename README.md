# BiofilmModel
## Code Execution
In order to run this model, open a case file in Matlab, input desired parameters, and run. Case file names are formatted as case_(name of case). If a new case file must be made, it is typically easiest to copy the format of an existing case file and update the parameters accordingly.

## Case Input Parameters
Open the case file and input all desired parameters.

param.instantaneousDiffusion – When true, uses biofilmdiffusion_fd to calculate instantaneous diffusion. Some cases run faster with this on.  
param.sourceTerm – Adds an extra source term to the growth equations. Was primarily used for dead cells that had a source term that depended on the death of live cells.  
param.Pulse – Turns on pulsed inlet for substrates.  
param.outPeriod – Sets how often the plots update. Useful for optimizing runtime.  
Param.Title – Adds title for plots.  
param.SNames – Defines substrate names for plot legends.  
param.XNames – Defines particulate names for plot legends.  
param.Xo – Particulate initial conditions in the tank. Formatted as a (# of particulates) x 1 array.  
param.So – Substrate initial conditions in the tank. Formatted as a (# of substrates) x 1 array.  
param.phibo – Particulates in the biofilm initial volume fractions. Must include values between 0 and 1 that sum to less than 1. Formatted as a (# of particulates) x 1 array.  
param.Sbo – Substrates in the biofilm initial conditions. Formatted as a (# of substrates) x 1 array.  
param.z – Grid points to solve over. Calculated as a linspace from 0 to the initial thickness with Nz points.  
param.dz – Space between the grid points. Calculated by taking the difference of two grid points.  
param.Daq – Substrate diffusion through the boundary layer. Formatted as a (# of substrates) x 1 array.  
param.De – Substrate diffusion through the biofilm. Formatted as a (# of substrates) x 1 array.  
param.rho – Particulate densities. Formatted as a (# of particulates) x 1 array.   
param.Yxs – Yield coefficient. Formatted as a (# of substrates) x (# of particulates).  
param.light – Equation for light intensity over time. Currently set as a truncated sin wave.  
param.mu – Growth rate equations for each individual particulate. Formatted as a (# of particulates) x 1 array of equations.  
param.phi_tot – Sum of particulate volume fractions. Calculated by taking the sum of phibo.  
param.Ns – Number of substrates. Found using the size of the So array.   
param.Nx – Number of particulates. Found using the size of the Xo array.  

The following function can be found at the bottom of the solver:  
function param = check_param(param)  

The purpose of this function is to check that all parameters have been entered correctly. If not, the function will either throw an error or enter filler parameters. For example, the neutralization parameter is set to false if it is not entered and yeild coefficients of 0 are changed to inf.


### Time Parameters
Determine the time parameters of the case being run.
The entire duration, the time-step interval size, and the tolerance for 
time-step convergence must all be defined.

tFin - duration [days]  
dt - time-step interval size [days]  
ttol - tolerance for time-step conversion  
outPeriod - frequency of plot updates


### Growth Rate Models
The growth rate model is determined by the equations entered into the param.mu variable. When updating growth model, param.Yxs will also have to be updated to reflect the new yeild coefficients.  


