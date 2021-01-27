% %% Week 1: Numerical Solution to First Order Incomplete ODE dx/dt
% clc; clear;
% 
% %Conditions
% v=100; %m^3
% mew=.00835; %growth coefficient
% Qp=1; %m^3/s %flow rate through tank
% t=linspace(1,2000,200); %seconds
% dt=t(2)-t(1); %delta time
% 
% %Preallocate
% x=zeros(1,length(t)); %Concentration of biomass in liquid
% N=length(x)-1; %length of for loop
% x(1)=.1; %Initial Condition
% 
% %Equation to solve ... dx/dt=mew*x(t)-(Q'*x(t))/V
% dxdt = @(x) x*(mew-(Qp/v))
% 
% %Eulers Method
% for i=1:N
%     x(i+1)=x(i)+dt*dxdt(x(i));
%     
% end
% 
% %Plot
% figure(1)
% plot(t,x)
% xlabel('time')
% ylabel('concentration of biomass in tank')

%% Week 2: Part of Solution to Second Order Diffusion (of Substrate Inside Biofilm)
clc;
clear;
% Future: Impliment time dependent concentrations


%Steady State

%Steady State Test Conditions
mewmax=20; %d-1	maximum specific growth rate
Km=3; %g m-3	Monod half saturation coefficient
Yxs=0.5; %gx gs-1	biomass yield coefficient on substrate
De=5.00E-05; %m2 d-1	effective diffusion coefficient of substrate in biofilm
Xb=20000; %g m-3	biomass density in biofilm
Lf=4.00E-04; %m	biofilm thickness
So=25; %g m-3	bulk fluid substrate concentration

%Setup for solution
N=100; %int64(Lf/dz); %necessary points (for full depth of biofilm)

%Grid of Biofilm Depth
z=linspace(0,Lf,N); %m
dz=z(2)-z(1); %m

%Substrate in Biofilm Boundary Conditions
S=zeros(1,N);
S(1)=0;
S(N)=So;

%Iterations
for iter=1:10000
    for i=2:N-1
        S(i)=(S(i+1)+S(i-1)-(mewmax*Xb*(dz^2))/(Yxs*De))/2; %Concentration of substrate at biofilm depth

        %Boundary Conditions After Iteration
        S(1)=0;
        S(end)=So;
        
        if S(i)<0
            S(i)=0;
        end
    end

end

figure(1); clf(1)
plot(z,S)
xlabel('Depth of Biofilm [m]','fontsize',20)
ylabel('Substrate Concentration within Biofilm [g/m^3]','fontsize',20)
