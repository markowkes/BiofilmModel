function [param]=case1()
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
tFin=5;   %[days]
dt  =1e-2; %time interval between calculations
ttol=1e-8; %tolerance for timestep conversion

%% Growth Rate Models
%Options
% model=1; %('Linear Growth Rate Equation')
% model=2; %('Monod Growth Rate Equation')
% model=3; %('Double Monod Growth Rate Equation')
% model=4; %('Inhibition Growth Rate Equation')
% model=5; %('None')
% model=[3 3]%('Double Monod Growth Rate for two Substrates

%% Biofilm
Nz=50; %Linear grid points to describe biofilm
dtol=1e-10; %tolerance for substrate diffusion convergence

%% Tank Geometry
L=0.5; %[m]
W=0.5; %[m]
H=0.5; %[m]

%% Frequency of Plots
outFreq=2000; %Number of steps between plot updates.

%% Substrate Dependant Properties
%model = [1 1; 1 1]; % j, k growth model for jth particulate based on kth substrate 

%% Constants
%mumax=[2000 2000];
%Km   =[2500 2; 2500 3]; %[3 3 3 3 3 3000 3 2500 3; 3 3 3 3 3 3000 3 2500 3];
%Km(:,:,2)=[3 1; 3 3];
Yxs  =[0.5];
V    =[0.1];
Q    =[1];
A    =(V/H)+2*((V/L)+(V/W));
Sin  =[25];
So   =[25];
Ns   =size(So, 1);
phio =[0.2];
xo   =[10];
Daq  =[4.0E-5];
De   =[1.0E-5];
Xb   =[20000];
rho  =[1.0E5];
Lfo  =[5.0E-6];
LL   =[1.00E-4];
Kdet =[1900];


% Growthrates for each biomass species 
% param.mu=@(j,S,param) [
%     (2000*S(1))./(2500)
%     (2000*S(2))./(2500)
%     ];
mu{1}=@(S,param) (20*S(1))./(3+S(1));

param.mu=mu;


%% Index variables under structure "param"
param.tFin   =tFin;
param.dtmax  =dt;
param.ttol   =ttol;

%param.model  =model;
param.Nx     =length(xo);

param.Nz     =Nz;
param.dtol   =dtol;

param.outFreq=outFreq;

%param.mumax  =mumax;
%param.Km     =Km;
param.Yxs    =Yxs(:);
param.V      =V;
param.Q      =Q;
param.A      =A;
param.Sin    =Sin(:);
param.So     =So(:);
param.Ns     =Ns;
param.phio   =phio;
param.phi_tot=sum(phio);
param.xo     =xo;
param.Daq    =Daq(:);
param.De     =De(:);
param.Xb     =Xb;
param.rho    =rho;
param.Lfo    =Lfo;
param.LL     =LL;
param.Kdet   =Kdet;
end 