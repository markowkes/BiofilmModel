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

