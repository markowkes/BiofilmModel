% Unit tests for biofilm model
function tests = TestCaseMultSubs
clear; clc
tests = functiontests(localfunctions);
end

%% Test when LL=0 
function test_diffusion_zeroLL(testCase)
% Run test
param=case1();
param.LL=0;
Nz=50;
Sbold=linspace(0,5,Nz);
S=10;
dz=1e-7;
t=0;
[Sb,~]=biofilmdiffusion_fd(Sbold,S,Nz,dz,t,param);
% Analyze result
figure(1); clf(1)
z=linspace(0,param.Lfo,Nz); %specify for plot
plot(z(end),Sb(end),'r*','Markersize',16)
hold on
plot(z,Sb,'black')
xl=xline(z(end),'--b','Biofilm Thickness','Fontsize',16);
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'center';
title('Substrate Concentration Profile')
xlabel('z')
ylabel('Sb(z)')
legend(sprintf('Sb = %3.3f [g/m^3]',Sb(end)),'Concentration Profile','location','Northwest','Fontsize',16)
set(gca,'Fontsize',20)

actSolution = Sb(end);
expSolution = S;
tol=1e-15;
verifyLessThan(testCase,abs(actSolution-expSolution),tol)
end

%% Test tank biomass concentration when no inflow Q
function test_tankenvironment_biomasssolution(testCase)
% Run Test
param=case1();
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
[~,~,x,~,~]=tankenvironment(t,x,S,Vdet,dt,bflux,param);
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
param=case1();
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
[~,~,~,S,~]=tankenvironment(t,x,S,Vdet,dt,bflux,param);
%Analyze Result
figure(1);clf(1);
plot(S)
actSolution=S;
expSolution=So;
tol=1e-1;
verifyLessThan(testCase,abs(actSolution-expSolution),tol)
end

%% Test Diffusion
function test_mult_subs_diffusion(testCase)
% Run Test
param.mumax=2000;
param.Km=2500;
param.Yxs=[0.5 -0.278];
S=[25; 2];
param.Daq=[4e-5 6e-5];
param.De =[1e-5 1.5e-5];
Lf=5e-6;
param.LL=0;
param.Xb=20000;
param.dtol=1e-12;
param.model=[1,1];
Yps     =1.8;
param.Ns=2;
param.Nx=1;
Xbo=1;

% Growthrates for each biomass species 
param.mu{1}=@(S,param) (param.mumax*S(1))./(param.Km);

figure(1); clf(1); hold on
figure(3); clf(3); hold on

Nzs=[10,25,50,100,200]; %Grid sizes to test
error=zeros(2,length(Nzs)); %Preallocate


for i=1:length(Nzs)
    Nz=Nzs(i);
    param.Nz=Nz; 
    
    % Define Xb
    Xb=zeros(1,Nz)+param.Xb;
    
    z=linspace(0,Lf,Nz); %[m] Grid of Biofilm Depth
    dz=z(2)-z(1); %[m]
    %Sbold=linspace(0,S(1),Nz);
    Sbold=zeros(2,Nz);
    t=0;
    [Sb,~]=biofilmdiffusion_fd(Sbold,S,Xb,dz,t,param);
    figure(1);
    plot(z,Sb(1,:),'DisplayName',['Gridsize: ',num2str(Nz)])
    hold on
    
    Stest=S(2)+(param.Daq(1)*Yps/param.Daq(2))*(S(1)-Sb(1,:));
    figure(3);
    plot(z,Stest)
    hold on
    plot(z,Sb(2,:),'--')
    legend('Analytic','Numerical')
    
    
    % Analyze Result
    phi = sqrt(param.mumax*param.Xb*Lf*Lf/...
        (param.De(1)*param.Km*param.Yxs(1)));
    Sb_ana = S(1)*cosh(phi*z/Lf)/cosh(phi);
    
    % Error
    error(1,i)=mean(abs(Sb(1,:)-Sb_ana));
    error(2,i)=mean(abs(Sb(2,:)-Stest));
    

end
figure(1)
plot(z,Sb_ana,'--','DisplayName','Analytic')
xlabel('z')
ylabel('Sb(z)')
title('Substrate Profiles within Biofilm')
legend('Location','Northwest')
set(gca,'Fontsize',20)

figure(2); clf(2)
loglog(Nzs,error(1,:),'-o')
hold on
loglog(Nzs,Nzs.^-1,'--')

xlabel('Number of grid points')
ylabel('Error')
set(gca,'Fontsize',20)

figure(3)
xlabel('z')
ylabel('Sb(z)')
title('Product Profiles within Biofilm')

% Pass/fail
tol=1e-2;
verifyLessThan(testCase,min(error),tol)
end

%% Test steady-state with large diffusivities
% Substrate concentration is driven to a constant in this case

function test_steadystate(testCase)
tc=7 ; %number of case, A corresponds to 1, B corresponds to 2....
param=case1(1); %structure variables are stored in
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
S=1; % Initial guess

tol=1e-12;
error=1;
while abs(error)>tol
    Lf=mu(S,param)/Kdet;
    Vdet=mu(S,param)*Lf;
    x=Yxs*(Sin-S);
    LHS=Q*x;
    RHS=Vdet*A*Xb+mu(S,param)*x*V;
    
    error=LHS-RHS;
    S=S+0.001*error;
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
param=case1();
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
maxError=max(abs(S-S_anal));  % Maximum Error
expTol=param.ttol;            % Expected Maximum Error
verifyLessThan(testCase,maxError,expTol)
end


%%

