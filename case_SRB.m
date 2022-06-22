% Multiple substrates 
clear; clc;

param.instantaneousDiffusion = true;
param.sourceTerm             = false;
param.Pulse                  = false;

% Time
param.tFin=30;   % Simulation time [days]
param.outPeriod=0.5; % Output period [days]

param.SNames = {'Oxygen', 'Sulfate', 'Sulfide'};
param.XNames = {'SRB'};

% Tank Geometry
param.V   = 0.01;        % Volume of tank
param.A   = 1;          % Surface area of biofilm
param.Q   = 1;          % Flowrate through tank
param.Xo  = 1;       % Tank particulate initial condition(s)
param.So  = [1; 48; 0]        % Tank substrate initial condition(s)

%param.SinPulse = [];
param.LL  = 2.00E-4;    % Boundary layer thickness

% Biofilm
param.Nz=25;          % Linear grid points to describe biofilm
param.phibo= 0.2;     % Biofilm particulates initial condition(s)
param.Sbo  = [0; 0; 0]     % Biofilm substrates initial condition(s)
param.Lfo  = 5.0E-6;  % Biofilm initial thickness
param.z  = linspace(0,param.Lfo,param.Nz);
param.dz = param.z(2) - param.z(1);

% Substance Constants
param.Daq  = [1.51E-4; 8.00E-5; 1.21E-4;]          % Substrate diffusion through boundary layer
param.De   = [6.8E-5;  4.00E-5; 6.04E-5;]        % Substrate diffusion through biofilm     
param.rho  = 2.5E5;            % Particulate densities
param.Kdet = 125;               % Particulates detachment coefficient
% Yield coeffficient
param.Yxs  = [0 0.584 -1.645]                % dX2/dS1  - Production of Oxygen

% Source term
param.b = 0.1;
X_Source{1}=@(S,X,param) 0;
param.X_Source=X_Source;

% Light term
param.I = 0.5;
param.diss = 1000;
param.Ylight = 2;
         
% Growthrates for each particulate
KmC2 = 20; KI = 1; mumaxC = 1.46;
param.light=@(t,z) heaviside(t-0.5)*max(0,param.I-(max(z)-z)*param.diss);
param.light=@(t,z) (cos(2*t)+1)*max(0,param.I-(max(z)-z)*param.diss); 
% mu{1}=@(S,X,t,z,param) (mumax*S(2,:))./(KmC2+S(2,:)).*(1./(1+S(1,:)/KI));
% param.mu=mu;
% param.light=light;

param.mu=@(S,X,t,z,param) [
    (mumaxC*S(2,:))./(KmC2+S(2,:)).*(1./(1+S(1,:)/KI));
];

% Computed parameters
param.phi_tot = sum(param.phibo);
param.Ns = size(param.So, 1);  % Number of substrates
param.Nx = size(param.Xo, 1);  % Number of substrates

param.Sin = [1; 48; 0];         % Substrates concentration(s) into tank
% Sin{1}.min   = 0;
% Sin{1}.max   = 50;
% Sin{1}.period= 15;
% Sin{1}.dur   = 10;
% % Sin{2}.max   = 0;
% 
% for k = 1:param.Ns
%     Sin{k}.f =@(theavi) (Sin{k}.max-Sin{k}.min)*(sum(heaviside(theavi)) ...
%     -sum(heaviside(theavi-Sin{k}.period+Sin{k}.dur)))+Sin{k}.min;
% end

% Tolerance
param.tol=1e-10;

%% Solver

% Call solver
[t,X,S,Pb,Sb,Lf]=solverBuiltIn(param);

% Plot solution
plotSolution(t,X,S,Pb,Sb,Lf,param)
