function MAINDRIVER
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
Lfo=4.00E-4; %biofilm thickness [m]
LL=Lfo/100; %thickness of boundary layer [m]
Lf_old=Lfo;
Co=So; %substrate concentration
Xb=20000; %g m^-3	biomass density in biofilm
De=5.00E-05; %m2 d^-1	effective diffusion coefficient of substrate in biofilm
Kdet=100/3600; % [1/ms] coefficient of detachment for biofilm

%Time Constraints
tFin=2; %[s]
dt=1e-3; %Interval
N=tFin/dt; %Number of steps

%Preallocation
t = zeros(1,N); %Time
x = zeros(1,N); %Biomass Concentration in bulk liquid
S = zeros(1,N); %Substrate in bulk liquid


%Initial Conditions
t(1)=0;
x(1)=xo;
S(1)=So;
Lf=Lfo;

for i = 1:N-1
   
    % Insert loop over spacial coordinate zeta for substrate diffusion,
    % particulates, biomass growth within biofilm, etc
    
    %Call on Biofilm Surface Substrate Concentration from 'Diffusion'
    [Sb,bflux,dz,z]=Diffusion(Lf,LL,S(i),mumax,Xb,Yxs,De);
    Cs=Sb(end);
    
    %Call on Biofilm Thickness and Vdet/Vg from 'BiofilmThickness_Fn'
    [Lf,Vdet,Vg]=BiofilmThickness_Fn(Sb,Lf_old,mu(Sb,mumax,Km),Kdet,mumax,Km,dt,dz);
    
    %Call on tank substrate/biomass concentration from 'tankenvironment'
    [S,x,t]=tankenvironment(xo,So,V,SA,Qdot,Sin,Vdet,mumax,Km,Yxs,Daq,LL,Cs,Co,Xb,dt,N);
    
    %Produce desired plots from outputs produced
    outputs(t,x,S,z,bflux,Sb);
    
   
    
    Lf_old=Lf;
    
end
end




