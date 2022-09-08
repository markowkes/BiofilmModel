% Multiple substrates 
clear; clc;

param.instantaneousDiffusion = true;
param.sourceTerm             = true;
param.Pulse                  = true;
param.neutralization         = true;

% Time
param.tFin=5;   % Simulation time [days]
param.outPeriod=0.1;

param.Title  = {'H2O2 Case'};
param.SNames = {'Substrate','H2O2'};
param.XNames = {'Live Particulate','Dead Particulate'};

% Tank Geometry
param.V   = 0.1;        % Volume of tank
param.A   = 1;          % Surface area of biofilm
param.Q   = 1;          % Flowrate through tank
param.Xo  = [1; 0];   % Tank particulate initial condition(s)
param.So  = [53; 0];   % Tank substrate initial condition(s)
param.LL  = 1.00E-5;    % Boundary layer thickness

% Biofilm
param.Nz   = 50;         %Linear grid points to describe biofilm
param.phibo= [0.08; 0]; % Biofilm particulates initial condition(s)
param.Sbo  = [25; 0];   % Biofilm substrates initial condition(s)
param.Lfo  = 5.0E-5;     % Biofilm initial thickness

% Substance Constants
param.Daq  =[5.70E-4; 1.73E-4];  % Substrate diffusion through boundary layer
param.De   =[2.57E-4; 7.79E-5];  % Substrate diffusion through biofilm     
param.rho  =[2.5E5; 2.5E5];    % Particulate densities
param.Kdet =1.5E4;             % Particulate detachment coefficient
% Yield coeffficients: new particulate/consumed substrate (can't be zero)
param.Yxs  =[0.5 0;        % [ dX1/dS1, dX1/dS2
             0   0];      %   dX2/dS1, dX2/dS2 ]

b = 0.5; param.kdis = [0; 2]; param.kB = [0; 5];
X_Source{1}=@(S,X,Pb,param) -param.kB(2)*param.kdis(2)*param.rho(1)*Pb(1,end);
X_Source{2}=@(S,X,Pb,param)  param.kB(2)*param.kdis(2)*param.rho(1)*Pb(1,end);
param.X_Source=X_Source;  

% Light term
param.I = 0.5;
param.diss = 0.05;
param.Ylight = 2;

KM = 2.55;
mumaxA = 7.5;
% Growthrates for each particulate
% mu{1}=@(S,X,t,z,param) (mumaxA*S(1,:))./(KM + S(1,:));
% mu{2}=@(S,X,t,z,param) 0*S(1,:);
% param.mu=mu;

param.mu=@(S,X,Lf,t,z,param) [
    (mumaxA*S(1,:))./(KM + S(1,:));
    0*S(1,:);
];


% Computed parameters
param.phi_tot = sum(param.phibo);
param.Ns = size(param.So, 1);  % Number of substrates
param.Nx = size(param.Xo, 1);  % Number of substrates

%param.Sin = [53; 0]; % Substrates concentration(s) into tank
% param.Sin{1}.min   = 0;
% param.Sin{1}.max   = 53;
% param.Sin{1}.period= 4;
% param.Sin{1}.dur   = 2;
% 
% param.Sin{2}.min   = 0;
% param.Sin{2}.max   = 53;
% param.Sin{2}.period= 4;
% param.Sin{2}.dur   = 2;

param.Sin.period = [4; 4];
param.Sin.min    = [0; 0];
param.Sin.max    = [53; 0];
param.Sin.dur    = [2; 2];

k = 1;
param.Sin.f{k} =@(theavi) 53;

k = 2;
param.Sin.f{k} =@(theavi) (param.Sin.max(k)-param.Sin.min(k))*(sum(heaviside(theavi)) ...
          -sum(heaviside(theavi-param.Sin.period(k)+param.Sin.dur(k))))+param.Sin.min(k);
%          -kB*kdis*param.rho(1)*Pb(1,end);


% Tolerance
param.tol=1e-2;

%% Solver

% Call solver
[t,X,S,Pb,Sb,Lf]=solverBuiltIn(param);

% Plot solution
plotSolution(t,X,S,Pb,Sb,Lf,param)
