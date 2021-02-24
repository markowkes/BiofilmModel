function MAINDRIVER
clear; clc

% Inputs

%Array Initial Conditions
xo=10; %[g/m^-3] initial biomass concentration in bulk liquid
So=25; %[g/m^-3] initial substrate concentration in bulk liquid

%Tank Parameters + Geometry
L=0.5; %[m]
W=0.5; %[m]
H=0.4; %[m]
V=L*W*H; %tank volume [m^3]
SA=(V/H)+2*((V/L)+(V/W)); %tank surface area [m^2] 
Qdot=1; %flow rate in/out [m^3]
Sin=25; %Inflow of substrates to tank, [g/m^3]

%Biofilm Parameters
mumax=20; %max specific growth rate
Km=3; %Monod half-saturation coefficient(growth transitions from sat. to linear)
Yxs=0.5; %ratio of substrate consumed to biomass produced
Daq=4.00e-5; %diffusion coefficient assumed at boundary [m/s^2]
Lfo=5.00E-6; %biofilm thickness [m]
Lf_old=Lfo;
LL=1.00e-4; %thickness of biofilm boundary layer [m]
Co=So; %substrate concentration
Xb=20000; %biomass density in biofilm [g/m^3]
De=1.00E-05; %effective diffusion coefficient of substrate in biofilm [m^2/d^1]
Kdet=1900; %coefficient of detachment for biofilm [m-1*s-1]

Nz=100; %Linear GridPoints in Biofilm
zo=linspace(0,Lfo,Nz); %[m] Grid of Biofilm Depth
dzo=zo(2)-zo(1); %[m]

%Initial Boundary Conditions (in Biofilm)
Sb=zeros(1,Nz);
Sb(end)=So; %initially assume boundary concentration = So

%Time Constraints
tFin=2; %[s]
dt=1e-3; %Interval
N=tFin/dt; %Number of steps
outFreq=10; %Number of steps between plot updates.

%Preallocation
t = zeros(1,N); %Time
x = zeros(1,N); %Biomass Concentration in bulk liquid
S = zeros(1,N); %Substrate in bulk liquid

bflux=zeros(1,N); %Boundary Layer Flux of Biofilm


%Initial Conditions
t(1)=0;
x(1)=xo;
S(1)=So;
Lf=Lfo;

% Initialize plots 
outIter=outFreq-1;
plots=0;

for i = 1:N-1
   
    % Insert loop over spacial coordinate zeta for substrate diffusion,
    % particulates, biomass growth within biofilm, etc
    
    %Call on Biofilm Surface Substrate Concentration from 'Diffusion'
    
    %Calculate Input Parameters
    Sbold=Sb; %Solve Diffusion with previous solution for efficiency
    z=linspace(0,Lf,Nz); %[m] Grid of Biofilm Depth
    dz=z(2)-z(1); %[m]
 
    [Cs,Sb,bflux(i)]=Diffusion(Sbold,Lf,LL,S(i),mumax,Xb,Yxs,De,Km,dz,Nz);
    
    %Call on Biofilm Thickness and Vdet/Vg from 'BiofilmThickness_Fn'
    [Lf,Vdet]=Lf(Sb,Lf_old,Kdet,mumax,Km,dt,dz);
    
    %Call on tank substrate/biomass concentration from 'tankenvironment'
    [t(i+1),x(i+1),S(i+1)]=tankenvironment(t(i),x(i),S(i),V,SA,Qdot,Sin,Vdet,mumax,Km,Yxs,Daq,LL,Cs,Co,Xb,dt);
    
    %Call on desired plots from 'outputs'
    outIter=outIter+1;
    if (outIter==outFreq)
        [plots] = outputs(t(1:i+1),x(1:i+1),S(1:i+1),z,bflux,Sb,plots);
        outIter=0;
    end
end
end




