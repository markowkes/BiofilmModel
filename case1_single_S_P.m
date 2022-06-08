% Multiple substrates 
clear; clc;

% Time
param.tFin=30;   % Simulation time [days]

% Tank Geometry
param.V   = 0.1;        % Volume of tank
param.A   = 1;          % Surface area of biofilm
param.Q   = 1;          % Flowrate through tank
param.Xo  = 10;         % Tank particulate initial condition(s)
param.So  = 100;         % Tank substrate initial condition(s)
param.Sin = 100;         % Substrates concentration(s) into tank
param.LL  = 1.00E-7;    % Boundary layer thickness

% Biofilm
param.Nz=50;            % Number of grid points to represent biofilm
param.phibo= 0.08;       % Biofilm particulates initial condition(s)
param.Sbo  = 0;         % Biofilm substrates initial condition(s)
param.Lfo  = 1.0E-5;    % Biofilm initial thickness

% Substance Constants
param.Yxs  = -2.646;       % Biomass yield coeffficient on substrate
param.Daq  = 4.0E-5;    % Substrate diffusion through boundary layer
param.De   = 6.9E-5;    % Substrate diffusion through biofilm     
param.rho  = 2.0E4;     % Particulate densities
param.Kdet = 20000;      % Particulate detachment coefficient

X_Source{1}=@(S,X,param) 0;
param.X_Source=X_Source;

% Growthrates for each particulate
mu{1}=@(S,X,t,z,param) (20*S(1,:))./(3+S(1,:));
param.mu=mu;  

% Computed parameters
param.phi_tot = sum(param.phibo);
param.Ns = size(param.So, 1);  % Number of substrates
param.Nx = size(param.Xo, 1);  % Number of substrates

% Tolerance
param.tol=1e-2;

%% Solver

% Call solver
[t,X,S,Pb,Sb,Lf]=solverBuiltIn(param);

% Plot solution
plotSolution(t,X,S,Pb,Sb,Lf,param)
