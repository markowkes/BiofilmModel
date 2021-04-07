function [param]=cases(tc)
%% This function outputs a set of constant parameters
% Each constant will have its own array of values for each case
% (A,B,C,D,E,F based on chronological order in the array).

%% Input Parameter Nomenclature
% tfin:
% dt:
% outFreq: 

% Nz:

% mumax: Maximum specific growth rate
% Km:    Monod Half saturation coefficient
% Yxs:   Biomass yield coeffficient on substrate
% V:     Volume of CSTR
% Q:     Flowrate
% SA:    Wetted Surface Area
% Sin:   Influent Substrate Concentration
% So:    Initial bulk fluid substrate concentration in tank
% xo:    Initial biomass concentration in tank
% Daq:   Diffusion coefficient of substrate in water
% De:    Effective diffusion coefficient of substrate in biofilm
% Xb:    Biomass density in biofilm
% Lfo:   Initial biofilm thickness
% LL:    Concentration boundary layer thickness
% Kdet:  Detachment rate coefficient
% SA:    Wetted Surface Area

%% Description of Test Cases
%Test Case G,7: elevated diffusion coefficients, so Sb profile approaches
%               constat/linear profile
%Test Case H,8: maximize growth rate, eliminate liquid layer LL to reach
%               constant concentration throughout biofilm Sb=S

% %% Test Case
% tc=7;

%% Time Constraints
tFin=30;   %[days]
dt  =1e-2; %time interval between calculations
ttol=1e-8; %tolerance for timestep conversion

%% Growth Rate Models
%Options
% model=1; %('Linear Growth Rate Equation')
model=2; %('Monod Growth Rate Equation')
% model=3; %('Double Monod Growth Rate Equation')
% model=4; %('Inhibition Growth Rate Equation')
% model=5; %('None')

%% Biofilm
Nz=50; %Linear grid points to describe biofilm
dtol=1e-10; %tolerance for substrate diffusion conversion

%% Tank Geometry
L=0.5; %[m]
W=0.5; %[m]
H=0.5; %[m]

%% Frequency of Plots
outFreq=2000; %Number of steps between plot updates.

%% Constants
mumax=[20 20 2 20 20 20 20 2000 0];
Km   =[3 3 3 3 3 3000 3 2500 3];
Yxs  =[0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5];
V    =[0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1];
Q    =[1 1 1 50 1 1 1 1 0.05];
SA   =(V(tc)/H)+2*((V(tc)/L)+(V(tc)/W));
Sin  =[25 25 25 25 25 25 25 25 25];
So   =[25 25 25 25 25 25 25 25 0];
xo   =[10 10 10 10 10 10 10 10 10 ];
Daq  =[4.0E-5 4.0E-5 4.0E-5 4.0E-5 4.0E-5 4.0E-5 4.0E2 4.00E-5 4.0E-5];
De   =[1.0E-5 1.0E-5 1.0E-5 1.0E-5 1.0E-5 1.0E-5 1.0E2 1.00E-5 1.0E-5];
Xb   =[20000 20000 20000 20000 20000 20000 20000 20000 20000];
Lfo  =[5.0E-6 3.0E-4 5.0E-6 5.0E-6 5.0E-6 5.0E-6 5.0E-6 5.00E-6 5.0E-6];
LL   =[1.0E-4 1.0E-4 1.0E-4 1.0E-4 1.0E-4 1.0E-4 1.0E-4 1.00E-7 1.0E-4];
Kdet =[1900 1900 1900 1900 190000 1900 1900 1900 1900];

%% Index variables under structure "param"
param.tFin   =tFin;
param.dt     =dt;
param.ttol   =ttol;

param.model  =model;

param.Nz     =Nz;
param.dtol   =dtol;

param.outFreq=outFreq;

param.mumax  =mumax(tc);
param.Km     =Km(tc);
param.Yxs    =Yxs(tc);
param.V      =V(tc);
param.Q      =Q(tc);
param.SA     =SA;
param.Sin    =Sin(tc);
param.So     =So(tc);
param.xo     =xo(tc);
param.Daq    =Daq(tc);
param.De     =De(tc);
param.Xb     =Xb(tc);
param.Lfo    =Lfo(tc);
param.LL     =LL(tc);
param.Kdet   =Kdet(tc);
end 