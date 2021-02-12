%% Heun's Method, Biomass and Substrate Concetration in Bulk Liquid 

clear; clc
% Inputs

%Array Initial Conditions
xo=10; %initial biomass concentration in bulk liquid
So=25; %initial substrate concentration in bulk liquid

%Tank Parameters + Geometry
L=0.5; %[m]
W=0.5; %[m]
H=0.4; %[m]
V=L*W*H; %tank volume [m^3]
SA=(V/H)+2*((V/L)+(V/W)); %tank surface area [m^2] 
Qdot=1; %[flow rate in/out, m^3]
Sin=25; %Inflow of substrates to tank, [g/m^3]

%Biofilm Parameters
mumax=20; %max specific growth rate
Km=3; %Monod half-saturation coefficient(growth transitions from sat. to linear)
Yxs=0.5; %ratio of substrate consumed to biomass produced
Daq=2e-5; %diffusion coefficient of water(assumed at boundary) [m/s^2]
Lf=4.00E-4; %biofilm thickness [m]
LL=Lf/100; %thickness of boundary layer [m]
Co=So; %substrate concentration
Xb=20000; %g m-3	biomass density in biofilm
De=5.00E-05; %m2 d-1	effective diffusion coefficient of substrate in biofilm

%Time Constraints
tFin=2; %[s]
dt=1e-3; %Interval
N=tFin/dt; %Number of steps

%S = @(x) Sin-(V/Qdot)*x+(V^2/Qdot); %Substrate concentration wrt biomass
%x = @(S) (Qdot*Sin/V)-(S*Qdot/V)+V; %Biomass wrt substrate concetration
%mu = @(x,S) ((mumax*S)/(Km+S)); %Monod Growth Kinetics

%Call on 'mu_function' for mu value
[mu] = mu_function(mumax,Km,S);

dsdt = @(x,t,S,Cs,mu) -((mu*x)/Yxs)+((Qdot*Sin)/V)-((Qdot*S)/V)-(SA*((Daq/LL)*(Co-Cs))); 
% ^^^Substrate Concentration Change wrt time, now also considers flux
% through boundary layer of biofilm
dxdt = @(x,t,S,Cs,mu) (mu-(Qdot/V))*x; %Biomass Concentration Change wrt time

%Preallocation
t = zeros(1,N); %Time
x = zeros(1,N); %Biomass Concentration in bulk liquid
S = zeros(1,N); %Substrate in bulk liquid

%Initial Conditions
t(1)=0;
x(1)=xo;
S(1)=So;

for i = 1:N-1
   
    % Insert loop over spacial coordinate zeta for substrate diffusion,
    % particulates, biomass growth within biofilm, etc
    
    %Calling Biofilm Surface Substrate Concentration from 'Diffusion'
    [Sb,bflux,flux]=Diffusion(Lf,S(i),mumax,Xb,Yxs,De);
    Cs=Sb(end);

    t(i+1) = t(i) + dt;
    
    xstar = x(i) + dt*dxdt(x(i),t(i),S(i),Cs,mu);
    Sstar = S(i) + dt*dsdt(x(i),t(i),S(i),Cs,mu);
    
    x(i+1) = x(i) + dt/2*(dxdt(x(i),t(i),S(i),Cs,mu)+dxdt(xstar,t(i+1),Sstar,Cs,mu));
    S(i+1) = S(i) + dt/2*(dsdt(x(i),t(i),S(i),Cs,mu)+dsdt(xstar,t(i+1),Sstar,Cs,mu)); 
    
end


%% plot
figure(2); clf(2)
plot(t,x)
hold on
plot(t,S)
title('Biomass and Substrate Concentrations For Filling/Draining Tank')
xlabel('Time')
ylabel('Amount of Biomass/Substrate in Tank')
legend('Biomass','Substrate')




