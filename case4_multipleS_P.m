% Multiple substrates 
clear; clc;

param.instantaneousDiffusion = true;
param.sourceTerm             = false;
param.Pulse                  = false;

param.Title  = {'Multiple S P Case'};
param.SNames = {'Substrate 1', 'Substrate 2'};
param.XNames = {'Bug 1', 'Bug 2'};

% Time
param.tFin=5;   % Simulation time [days]
param.outPeriod=1e-1;

% Tank Geometry
param.V   = 0.1;        % Volume of tank
param.A   = 1;          % Surface area of biofilm
param.Q   = 1;          % Flowrate through tank
param.Xo  = [10; 10];   % Tank particulate initial condition(s)
param.So  = [25; 25];   % Tank substrate initial condition(s)
param.LL  = 1.00E-4;    % Boundary layer thickness

% Biofilm
param.Nz   = 50;         %Linear grid points to describe biofilm
param.phibo= [0.2; 0.2]; % Biofilm particulates initial condition(s)
param.Sbo  = [25; 25];   % Biofilm substrates initial condition(s)
param.Lfo  = 5.0E-6;     % Biofilm initial thickness

% Substance Constants
param.Daq  =[4.0E-5; 6.0E-5];  % Substrate diffusion through boundary layer
param.De   =[1.0E-5; 1.5E-5];  % Substrate diffusion through biofilm     
param.rho  =[3.0E5; 3.0E5];    % Particulate densities
param.Kdet = 1900;             % Particulate detachment coefficient
% Yield coeffficients: new particulate/consumed substrate (can't be zero)
param.Yxs  =[0.5,   inf        % [ dX1/dS1, dX1/dS2
             inf, 0.278];      %   dX2/dS1, dX2/dS2 ]

b = 0;
X_Source{1}=@(S,X,Pb,param) -b*X(1,:);
X_Source{2}=@(S,X,Pb,param)  b*X(1,:);
param.X_Source=X_Source;  

% Light term
param.I = 0.5;
param.diss = 0.05;
param.Ylight = 2;
         
% Growthrates for each particulate
mumax = 2000; Km = 2500;
% mu{1}=@(S,X,param) (2000*S(1,:))./(2500);
% mu{2}=@(S,X,param) (2000*S(2,:))./(2500);
% param.mu=mu;

param.mu=@(S,X,Lf,t,z,param) [
    (mumax*S(1,:))./(Km);
    (mumax*S(2,:))./(Km);
];

% Computed parameters
param.phi_tot = sum(param.phibo);
param.Ns = size(param.So, 1);  % Number of substrates
param.Nx = size(param.Xo, 1);  % Number of substrates

param.Sin.f{1}= @(theavi) 25;
param.Sin.f{2}= @(theavi) 25;
param.Sin.period = [0; 0];

% Tolerance
param.tol=1e-2;

%% Solver

% Call solver
[t,X,S,Pb,Sb,Lf]=solverBuiltIn(param);

% Plot solution
plotSolution(t,X,S,Pb,Sb,Lf,param)

