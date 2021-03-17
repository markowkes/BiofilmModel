% Unit tests for biofilm model
function tests = biofilmTest
clear; clc
tests = functiontests(localfunctions);
end

%% Test flux matching
function test_diffusion_flux(testCase)
% Run test
param=cases(1);
Nz=50;
Sbold=linspace(0,5,Nz);
S=10;
dz=1e-7;
[~,Sb,bflux,flux]=biofilmdiffusion_fd(Sbold,S,Nz,dz,param);
% Analyze result
figure(1); clf(1)
plot(Sb)
actSolution = bflux-flux;
expSolution = 0;
tol=1e-12;
verifyLessThan(testCase,abs(actSolution-expSolution),tol)
end

%% Test when LL=0 
function test_diffusion_zeroLL(testCase)
% Run test
param=cases(1);
param.LL=0;
Nz=50;
Sbold=linspace(0,5,Nz);
S=10;
dz=1e-7;
[~,Sb,~,~]=biofilmdiffusion_fd(Sbold,S,Nz,dz,param);
% Analyze result
figure(1); clf(1)
plot(Sb)
actSolution = Sb(end);
expSolution = S;
tol=1e-15;
verifyLessThan(testCase,abs(actSolution-expSolution),tol)
end

%% Test solution of diffusion problem
function test_diffusion_solution(testCase)
% Run test
param=cases(1);
Nz=50;
Sbold=linspace(0,5,Nz);
S=10;
dz=1e-7;
[~,Sb,~,~]=biofilmdiffusion_fd(Sbold,S,Nz,dz,param);
% Analyze result
figure(1); clf(1)
plot(Sb)
actSolution = Sb(1);
expSolution = 3.921;
tol=1e-3;
verifyLessThan(testCase,abs(actSolution-expSolution),tol)
end

%% Test tank biomass concentration when no inflow Q
function test_tankenvironment_biomasssolution(testCase)
% Run Test
param=cases(1);
param.Q=0;
L=0.5; %[m]
W=0.5; %[m]
H=0.4; %[m]
SA=(param.V/H)+2*((param.V/L)+(param.V/W)); %tank surface area [m^2] 
tFin=20; %[days]
dt=1e-2; %Interval
N=tFin/dt; %Number of steps
Vdet=2.7244e-5;
Cs=0.7079;
Co=25;
xo=param.xo;
t = zeros(1,N); %Time
x = zeros(1,N); %Biomass Concentration in bulk liquid
S = zeros(1,N); %Substrate in bulk liquid
[~,x,~,~]=tankenvironment(t,x,S,SA,Vdet,dt,Cs,Co,param);
% Analyze result
figure(1); clf(1)
plot(x)
actSolution=x;
expSolution=0;
tol=1e-1;
verifyLessThan(testCase,abs(actSolution-expSolution),tol)
end

%% Test tank substrate concentration when no inflow Sin
function test_tankenvironment_substratesolution(testCase)
%Run Test
param=cases(1);
param.Sin=0;
L=0.5; %[m]
W=0.5; %[m]
H=0.4; %[m]
SA=(param.V/H)+2*((param.V/L)+(param.V/W)); %tank surface area [m^2] 
tFin=20; %[days]
dt=1e-2; %Interval
N=tFin/dt; %Number of steps
Vdet=2.7244e-5;
Cs=0.7079;
Co=25;
So=param.So;
t = zeros(1,N); %Time
x = zeros(1,N); %Biomass Concentration in bulk liquid
S = zeros(1,N); %Substrate in bulk liquid
[~,~,S,~]=tankenvironment(t,x,S,SA,Vdet,dt,Cs,Co,param);
%Analyze Result
figure(1);clf(1);
plot(S)
actSolution=S;
expSolution=0;
tol=1e-1;
verifyLessThan(testCase,abs(actSolution-expSolution),tol)
end
