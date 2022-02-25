% Multiple substrates 
clear; clc;

% Time
param.tFin=5;   % Simulation time [days]

% Tank Geometry
param.V   = 0.1;        % Volume of tank
param.A   = 1;          % Surface area of biofilm
param.Q   = 1;          % Flowrate through tank
param.Xo  = [10; 10];   % Tank particulate initial condition(s)
param.So  = [25; 25];   % Tank substrate initial condition(s)
param.Sin = [25; 25];   % Substrates concentration(s) into tank
param.LL  = 1.00E-4;    % Boundary layer thickness

% Biofilm
param.Nz=50; %Linear grid points to describe biofilm
param.phibo= [0.2; 0.2]; % Biofilm particulates initial condition(s)
param.Sbo  = [25; 25];   % Biofilm substrates initial condition(s)
param.Lfo  = 5.0E-6;    % Biofilm initial thickness

% Substance Constants
param.Yxs  =[0.5,   inf        % X1 yield coeffficient on substrates
             inf,  0.278];     % X2 yield coeffficient on substrates
param.Daq  =[4.0E-5; 6.0E-5];  % Substrate diffusion through boundary layer
param.De   =[1.0E-5; 1.5E-5];  % Substrate diffusion through biofilm     
param.rho  =[3.0E5; 3.0E5];    % Particulate densities
param.Kdet = 1900;             % Particulate detachment coefficient
% Yield coeffficients: new particulate/consumed substrate (can't be zero)
param.Yxs  =[0.5,   inf        % [ dX1/dS1, dX1/dS2
             inf, 0.278];      %   dX2/dS1, dX2/dS2 ]

% Computed parameters
param.phi_tot = sum(param.phibo);
param.Ns = size(param.So, 1);  % Number of substrates
param.Nx = size(param.Xo, 1);  % Number of substrates

% Growthrates for each particulate
mu{1}=@(S,param) (2000*S(1))./(2500);
mu{2}=@(S,param) (2000*S(2))./(2500);
param.mu=mu;  

%% Solver

% Call solver
[t,X,S,Pb,Sb,Lf]=solverBuiltIn(param);

% Plot solution
plotSolution(t,X,S,Pb,Sb,Lf,param)

