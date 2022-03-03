% Unit tests for biofilm model
function tests = biofilmTest
clear; clc
tests = functiontests(localfunctions);
end


%% Test when LL=0 
function test_diffusion_zeroLL(testCase)
% Run test
param.Daq  = 4.0E-5;    % Substrate diffusion through boundary layer
param.De   = 1.0E-3;    % Substrate diffusion through biofilm     
param.Yxs  = 0.5;       % Biomass yield coeffficient on substrate
param.LL=1e-2;
param.Ns=1;
param.Nx=1;
param.Nz=50;
param.tol=1e-2;
param.Lfo=50e-6;
mu{1}=@(S,param) (20*S(1,:))./(3+S(1,:));
param.mu=mu;  

S=10;
Xb=zeros(1,param.Nz)+20000;
param.dtol=1e-2;

grid.z=linspace(0,param.Lfo,param.Nz); %specify for plot
grid.dz = grid.z(2)-grid.z(1);
% Old method
[Sbold,~]=biofilmdiffusion_fd_old(S,Xb,param,grid);
figure(1); clf(1)
plot(grid.z(end),Sbold(end),'r*','Markersize',16)
hold on
plot(grid.z,Sbold,'r--')

% New method
grid.z=linspace(0,param.Lfo,param.Nz+1); %specify for plot
grid.zm=0.5*(grid.z(1:end-1)+grid.z(2:end));
grid.dz = grid.z(2)-grid.z(1);
[Sb,~]=biofilmdiffusion_fd(S,Xb,param,grid);
Stop = ((param.Daq*S)*grid.dz + (2*param.De*Sb(:,end))*param.LL) ...
    /(param.Daq*grid.dz + (2*param.De)*param.LL);
plot(grid.zm,Sb,'black')
plot(grid.z(end),Stop,'bo','Markersize',16)

% Annotate
xl=xline(grid.z(end),'--b','Biofilm Thickness','Fontsize',16);
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'center';
title('Substrate Concentration Profile')
xlabel('z')
ylabel('Sb(z)')
legend(sprintf('Sb = %3.3f [g/m^3]',Sb(end)),'Concentration Profile','location','Northwest','Fontsize',16)
set(gca,'Fontsize',20)

actSolution = Stop;
expSolution = S;
error=abs(actSolution-expSolution);
tol=1e-15;
verifyLessThan(testCase,error,tol)
end

%% Test tank biomass concentration when no inflow Q
function test_tankenvironment_biomasssolution(testCase)
% Run Test
param=cases(1);
param.Q=0;
tFin=20; 
dt=1e-2; 
N=tFin/dt; 
t=0; 
xo=param.xo;
x=param.xo;
S=param.So;
bflux=0;
Vdet=0;
[~,x,~,~]=tankenvironment(t,x,S,Vdet,dt,bflux,param);
% Analyze result
figure(1); clf(1)
plot(x)
actSolution=x;
expSolution=xo;
tol=1e-1;
verifyLessThan(testCase,abs(actSolution-expSolution),tol)
end

%% Test tank substrate concentration when no inflow Q
function test_tankenvironment_substratesolution(testCase)
%Run Test
param=cases(1);
param.Q=0;
tFin=20; 
dt=1e-2; 
N=tFin/dt; 
So=param.So;
t=0; 
x=param.xo;
S=param.So;
bflux=0;
Vdet=0;
[~,~,~,S,~]=tankenvironment(t,x,S,Vdet,Xb,dt,bflux,param);
%Analyze Result
figure(1);clf(1);
plot(S(1))
actSolution=S(1);
expSolution=So(1);
tol=1e-1;
verifyLessThan(testCase,abs(actSolution-expSolution),tol)
end

%% Test diffusion
% Note that the analytic solution assumes mu=mumax*S/Km
% Method shows first order convergence
function test_diffusion(testCase)
% Run Test
param.mumax=2000;
param.Km=2500;
param.Yxs=0.5;
S=25;
param.Daq=4e-5;
param.De =1e-5;
Lf=5e-6;
param.LL=0;
param.Xb=20000;
param.dtol=1e-12;
param.model=1;
param.Ns=1;
param.Nx=1;

% Growthrates for each biomass species 
param.mu{1}=@(S,param) (param.mumax*S(1))./(param.Km);

figure(1); clf(1); hold on
Nzs=[10,50,100,1000,2000]; %Grid sizes to test
error=zeros(1,length(Nzs)); %Preallocate

for i=1:length(Nzs)
    Nz=Nzs(i);
    param.Nz=Nz;
    
    % Define Xb
    Xb=zeros(1,Nz)+param.Xb;
    
    z=linspace(0,Lf,Nz); %[m] Grid of Biofilm Depth
    dz=z(2)-z(1); %[m]
    Sbold=linspace(0,S,Nz);
    t=0;
    [Sb,~]=biofilmdiffusion_fd(Sbold,S,Xb,dz,t,param);
    plot(z,Sb)
    
    % Analyze Result
    phi = sqrt(param.mumax*param.Xb*Lf*Lf/...
        (param.De*param.Km*param.Yxs));
    Sb_ana = S*cosh(phi*z/Lf)/cosh(phi);
    
    % Error
    error(i)=mean(abs(Sb-Sb_ana));
end
plot(z,Sb_ana,'--')
xlabel('z')
ylabel('Sb(z)')
title('Substrate Profiles within Biofilm')
legend(sprintf('Gridsize:%5.0f',Nzs(1)),...
       sprintf('Gridsize:%5.0f',Nzs(2)),...
       sprintf('Gridsize:%5.0f',Nzs(3)),...
       sprintf('Gridsize:%5.0f',Nzs(4)),...
       sprintf('Gridsize:%5.0f',Nzs(5)),...
       'Analytic','Location','Northwest')
set(gca,'Fontsize',20)

figure(2); clf(2)
loglog(Nzs,error,'-o')
hold on
loglog(Nzs,Nzs.^-1,'--')

xlabel('Number of grid points')
ylabel('Error')
set(gca,'Fontsize',20)

% Pass/fail
tol=1e-2;
verifyLessThan(testCase,min(error),tol)
end

%% Test steady-state with large diffusivities
% Substrate concentration is driven to a constant in this case

function test_steadystate(testCase)
tc=7 ; %number of case, A corresponds to 1, B corresponds to 2....
param=cases(tc); %structure variables are stored in
% Run simulation
[~,xsim,Ssim,Lfsim]=MAINDRIVER(param);

% Analytic solution
Yxs=param.Yxs;
Sin=param.Sin;
Kdet=param.Kdet;
V=param.V;
A=param.A;
Q=param.Q;
Xb=param.Xb;
S=[0,0]; % Initial guess

tol=1e-12;
error=1;
while abs(error)>tol
    Lf=mu(1,S(1,:),param)/Kdet;
    Vdet=mu(1,S(1,:),param)*Lf;
    x=Yxs(1)*(Sin(1,:)-S(1,:));
    LHS=Q*x;
    RHS=Vdet*A*Xb+mu(1,S(1,:),param)*x*V;
    
    error=LHS-RHS;
    S(1,:)=S(1,:)+0.001*error(1);
end
% Compare solution
fprintf('S      =%16.12f, %16.12f g/m^3 \n',S,Ssim(end))
fprintf('x      =%16.12f, %16.12f g/m^3 \n',x,xsim(end))
fprintf('Lf     =%16.12f, %16.12f mu m \n',Lf*1e6,Lfsim(end)*1e6)

% Pass/fail
tol=.01; 
verifyLessThan(testCase,max((S-Ssim(end))/S,(x-xsim(end)))/x,tol)
end

%% Test time dynamic of tank environment calculations for dt
function test_timedynamicsdt(testCase)
% run test

% Simulation with no biomass (mu=0)
param=cases(9);
[t,~,S,~]=MAINDRIVER(param);

% Analytic Solution
param=cases(9);
Q=param.Q;
V=param.V;
Sin=param.Sin;
So=param.So;
S_anal=Sin*(1-exp(-Q/V*t))+So;

% Compare Simulation and Analytic
figure(1); clf(1)
plot(t,S)
hold on
plot(t,S_anal,'--')
legend('Simulation','Analytic')
title('Convergence of Simulated vs Analytical Methods')
ylabel('Output')
xlabel('Time Iteration')

%Analyze Result
maxError=max(abs(S(1,:)-S_anal(1,:)));  % Maximum Error
expTol=param.ttol;            % Expected Maximum Error
verifyLessThan(testCase,maxError,expTol)
end
