% Unit tests for biofilm model
function tests = biofilmTest
clear; clc
tests = functiontests(localfunctions);
end

%% Test flux matching
function test_diffusion(testCase)
% Run test
param=cases(1);
Nz=50;
Sbold=linspace(0,5,Nz);
S=10;
dz=1e-3;
[~,~,bflux,flux]=biofilmdiffusion_fd(Sbold,S,Nz,dz,param);
% Analyze result
actSolution = bflux-flux;
expSolution = 0;
tol=1e-15;
verifyLessThan(testCase,abs(actSolution-expSolution),tol)
end

%% Test flux when LL=0 
function test_zeroLL(testCase)
% Run test
param=cases(1);
param.LL=0;
Nz=50;
Sbold=linspace(0,5,Nz);
S=10;
dz=1e-3;
[~,Sb,~,~]=biofilmdiffusion_fd(Sbold,S,Nz,dz,param);
% Analyze result
actSolution = max(abs(Sb-S));
expSolution = 0;
tol=1e-15;
verifyLessThan(testCase,abs(actSolution-expSolution),tol)
end