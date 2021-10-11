% Unit tests for biofilm model
function tests = TestCaseMultSubs
clear; clc
tests = functiontests(localfunctions);
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
param.model=1;
Yps     =1.8;

figure(1); clf(1); hold on
Nzs=[10,50,100,1000,2000]; %Grid sizes to test
error=zeros(2,length(Nzs)); %Preallocate

for i=1:length(Nzs)
    Nz=Nzs(i);
    z=linspace(0,Lf,Nz); %[m] Grid of Biofilm Depth
    dz=z(2)-z(1); %[m]
    %Sbold=linspace(0,S(1),Nz);
    Sbold=zeros(2,Nz);
    t=0;
    [Sb,~]=biofilmdiffusion_fd(Sbold,S,Nz,dz,t,param);
    figure(1)
    plot(z,Sb(1,:))
    
    Stest=S(2)+(param.Daq(1)*Yps/param.Daq(2))*(S(1)-Sb(1,:));
    figure(3)
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

