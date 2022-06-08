% Multiple substrates 
clear; clc;

% Time
param.tFin=30;   % Simulation time [days]

param.SNames = {'Oxygen', 'Sulfate', 'Sulfide'};
param.XNames = {'SOB', 'SRB'};

% Tank Geometry
param.V   = 0.01;        % Volume of tank
param.A   = 1;          % Surface area of biofilm
param.Q   = 10;          % Flowrate through tank
param.Xo  = [1; 1]       % Tank particulate initial condition(s)
param.So  = [8.6; 48; 0]        % Tank substrate initial condition(s)

%param.SinPulse = [];
param.LL  = 2.00E-4;    % Boundary layer thickness

% Biofilm
param.Nz=25;          % Linear grid points to describe biofilm
param.phibo= [0.1; 0.1]    % Biofilm particulates initial condition(s)
param.Sbo  = [0; 0; 0]     % Biofilm substrates initial condition(s)
param.Lfo  = 5.0E-6;  % Biofilm initial thickness
param.z  = linspace(0,param.Lfo,param.Nz);
param.dz = param.z(2) - param.z(1);

% Substance Constants
param.Daq  = [1.51E-4; 8.00E-5; 1.21E-4;]          % Substrate diffusion through boundary layer
param.De   = [6.8E-5;  4.00E-5; 6.04E-5;]        % Substrate diffusion through biofilm     
param.rho  = [2.5E5;  2.5E5;]            % Particulate densities
param.Kdet = 120;               % Particulates detachment coefficient
% Yield coeffficient
param.Yxs  = [0.0581 0 0.09; 0 0.584 -1.645]                % dX2/dS1  - Production of Oxygen

% Source term
param.b = 0.1;
X_Source{1}=@(S,X,param) 0;
X_Source{2}=@(S,X,param) 0;
param.X_Source=X_Source;

% Light term
param.I = 0.5;
param.diss = 1000;
param.Ylight = 2;
         
% Growthrates for each particulate
KmB1 = 1; KmB3 = 11;  KmC2 = 20;  KI = 10;
mumaxB = 6.72;  mumaxC = 1.46;
light=@(t,z) (cos(2*t)+1)*max(0,param.I-(max(z)-z)*param.diss);
mu{1}=@(S,X,t,z,param) (mumaxB*S(1,:))./(KmB1+S(1,:)).*(S(3,:))./(KmB3+S(3,:));
mu{2}=@(S,X,t,z,param) (mumaxC*S(2,:))./(KmC2+S(2,:)).*(1./(1+S(1,:)/KI));
param.mu=mu;
param.light=light;

% Computed parameters
param.phi_tot = sum(param.phibo);
param.Ns = size(param.So, 1);  % Number of substrates
param.Nx = size(param.Xo, 1);  % Number of substrates

param.Sin = [8.6; 48; 0];         % Substrates concentration(s) into tank
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
