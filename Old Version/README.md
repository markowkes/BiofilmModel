# BiofilmModel
## Code Execution
In order to run this model, open the MAINDRIVER.m file in Matlab and type 
”MAINDRIVER(number)” in the command window. The number within the 
parentheses will correspond to the desired case to be run, which can be 
inputted in the cases.m file.

## Case Input Parameters
Open the cases.m file and input all given parameters into defined arrays 
under the "Constants" section. The array position for the particular case 
that is being run will be used as the only input when MAINDRIVER.m is run.

mumax: Maximum specific growth rate [1/days] 
Km:    Monod Half saturation coefficient [g/m^3]
Kib:   Inhibition Coefficient [g/m^3]
Yxs:   Biomass yield coeffficient on substrate [g g/s]
V:     Volume of CSTR [m^3]
Q:     Flowrate [m^3/day]
A:     Wetted Surface Area [m^2]
Sin:   Influent Substrate Concentration [g/m^3]
So:    Initial bulk fluid substrate concentration in tank [g/m^3]
xo:    Initial biomass concentration in tank [g/m^3]
Xb:    Biomass density in biofilm [g/m^3]
Daq:   Diffusion coefficient of substrate in water [m^2/day]
De:    Effective diffusion coefficient of substrate in biofilm [m^2/day]
Lfo:   Initial biofilm thickness [m]
LL:    Concentration boundary layer thickness [m]
Kdet:  Detachment rate coefficient [1/(m days)]

### Time Parameters
The next step is to determine the time parameters of the case being run.
The entire duration, the time-step interval size, and the tolerance for 
time-step convergence must all be defined.

tFin:  duration [days]
dt:    time-step interval size [days]
ttol:  tolerance for time-step conversion

### Growth Rate Models
The next step is to determine the desired growth rate model for this case.
This can be done in the "Growth Rate Models" section where there are five 
available models.

model 1: Linear Growth Rate
model 2: Monod Growth Rate 
model 3: Double Monod Growth Rate
model 4: Inhibition Growth Rate
model 5: None

### Biofilm Gridsize
Now the biofilm gridsize must be defined. The gridsize will determine the
precision of the solution this model produces, as a higher gridsize will 
produce a higher precision in the solutions. The tolerance for diffusion 
convergence must also be defined.

Nz:   Linear grid points to describe biofilm
dtol: tolerance for substrate diffusion convergence

### Tank Geometry
This step requires the size and geometry of the tank for the controlled 
stirred tank reactor (CSTR) to be defined.

L:    Length [m]
W:    Width [m]
H:    Height [m]

### Frequency of Plots
This section requires determination for the amount of steps between
each update in the plots being produced. When high volumes of iterations 
are being used, it makes sense to increase the time between plot updates 
to inmprove the efficiency of the model.

outfreq: Number of steps between plot updates
