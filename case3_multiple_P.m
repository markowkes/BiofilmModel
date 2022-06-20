% Multiple substrates 
clear; clc;

param.instantaneousDiffusion = true;
param.sourceTerm             = false;
param.Pulse                  = false;

param.SNames = {'Substrate'};
param.XNames = {'Bug 1', 'Bug 2'};

% Time
param.tFin=100;   % Simulation time [days]
param.outPeriod=1e-1;

% Tank Geometry
param.V   = 0.1;        % Volume of tank
param.A   = 1;          % Surface area of biofilm
param.Q   = 1;          % Flowrate through tank
param.Xo  = [10; 0];    % Tank particulate initial condition(s)
param.So  = 25;         % Tank substrate initial condition(s)

param.SinPulse = 30;
param.LL  = 1.00E-5;    % Boundary layer thickness

% Biofilm
param.Nz=12; %Linear grid points to describe biofilm
param.phibo= [0.08;0]; % Biofilm particulates initial condition(s)
param.Sbo  = 0;         % Biofilm substrates initial condition(s)
param.Lfo  = 5.0E-6;    % Biofilm initial thickness

% Substance Constants
param.Daq  = 1.38E-4;          % Substrate diffusion through boundary layer
param.De   = 6.9E-5;          % Substrate diffusion through biofilm     
param.rho  =[2.5E5; 2.5E5];   % Particulate densities
param.Kdet = 1980;            % Particulates detachment coefficient
% Yield coeffficient
param.Yxs  =[0.378                % dX1/dS1
             0];              % dX2/dS1

% Source term
param.b = 0.1;
X_Source{1}=@(S,X,param) -param.b*X(1,:);
X_Source{2}=@(S,X,param)  param.b*X(1,:);
param.X_Source=X_Source;

% Light term
param.I = 0.5;
param.diss = 50000;
param.Ylight = 2;
         
% Growthrates for each particulate
Km = 1; mumax = 2;
% mu{1}=@(S,X,t,z,param) (mumax*S(1,:))./(Km+S(1,:));
% mu{2}=@(S,X,t,z,param)  zeros(size(S(1,:)));
% param.mu=mu;  

param.mu=@(S,X,t,z,param) [
    (mumax*S(1,:))./(Km+S(1,:));
    zeros(size(S(1,:)));
];

% Computed parameters
param.phi_tot = sum(param.phibo);
param.Ns = size(param.So, 1);  % Number of substrates
param.Nx = size(param.Xo, 1);  % Number of substrates

param.Sin = [25];         % Substrates concentration(s) into tank
Sin{1}.min   = 0;
Sin{1}.max   = 50;
Sin{1}.period= inf;
Sin{1}.dur   = 10;
% Sin{2}.max   = 0;

for k = 1:param.Ns
    Sin{k}.f =@(theavi) (Sin{k}.max-Sin{k}.min)*(sum(heaviside(theavi)) ...
    -sum(heaviside(theavi-Sin{k}.period+Sin{k}.dur)))+Sin{k}.min;
end

% Tolerance
param.tol=1e-10;

%% Solver

% Call solver
[t,X,S,Pb,Sb,Lf]=solverBuiltIn(param);

% Plot solution
plotSolution(t,X,S,Pb,Sb,Lf,param)

