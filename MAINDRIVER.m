 %MAINDRIVER
clear; clc
tic
%% Rewriting inputs section to run through test cases
% TC=6;
num=1 ; %number of case, a corresponds to 1, b corresponds to 2....
param=cases(num); %all variables under structure "param"

% for j=1:TC %Number of Test Case
%     mumax=inputvariables(1,j);
%     Km   =inputvariables(2,j);
%     Yxs  =inputvariables(3,j);
%     V    =inputvariables(4,j);
%     Q    =inputvariables(5,j);
%     A    =inputvariables(6,j);
%     Sin  =inputvariables(7,j);
%     So   =inputvariables(8,j);
%     xo   =inputvariables(9,j);
%     Daq  =inputvariables(10,j);
%     De   =inputvariables(11,j);
%     Xb   =inputvariables(12,j);
%     Lf  =inputvariables(13,j);
%     LL   =inputvariables(14,j);
%     Kdet =inputvariables(15,j);

% %% Inputs
% 
% %Array Initial Conditions
% xo=10; %[g/m^-3] initial biomass concentration in bulk liquid
% So=25; %[g/m^-3] initial substrate concentration in bulk liquid

% Tank Parameters + Geometry
 L=0.5; %[m]
 W=0.5; %[m]
 H=0.4; %[m]
% V=L*W*H; %tank volume [m^3]
SA=(param.V/H)+2*((param.V/L)+(param.V/W)); %tank surface area [m^2] 
 
% Q=1; %flow rate in/out [m^3]
% Sin=25; %Inflow of substrates to tank, [g/m^3] 
% % Biofilm Parameters
% mumax=20; %max specific growth rate
% Km=3; %Monod half-saturation coefficient(growth transitions from sat. to linear)
% Yxs=0.5; %ratio of substrate consumed to biomass produced
% Daq=4.00e-5; %diffusion coefficient assumed at boundary [m/s^2]
% Xb=20000; %biomass density in biofilm [g/m^3]
% De=1.00e-5; %effective diffusion coefficient of substrate in biofilm [m^2/d^1]
% Kdet=1900; %coefficient of detachment for biofilm [m-1*s-1]
% Lf=5.00E-6; %biofilm thickness [m] 
% LL=1.00e-4; %thickness of biofilm boundary layer [m]

 Co=param.So; %substrate concentration

Nz=50; %Linear GridPoints in Biofilm
z=linspace(0,param.Lf,Nz); %[m] Grid of Biofilm Depth
dz=z(2)-z(1); %[m]

%Initial Boundary Conditions (in Biofilm)
Sb=zeros(1,Nz);
Sb(end)=param.So; %initially assume boundary concentration = So

% Time Constraints
tFin=2; %[s]
dt=1e-3; %Interval
N=tFin/dt; %Number of steps
outFreq=10; %Number of steps between plot updates.

% Preallocation
t = zeros(1,N); %Time
x = zeros(1,N); %Biomass Concentration in bulk liquid
S = zeros(1,N); %Substrate in bulk liquid
bflux=zeros(1,N); %Boundary Layer Flux of Biofilm Preallocate
flux=zeros(1,N); %Right hand side of power point equation to ensure matching flux

% Initial Conditions
t(1)=0;
x(1)=param.xo;
S(1)=param.So;

% Initialize plots 
outIter=outFreq-1;
plots=0;

%% Time Loop
for i = 1:N-1
    % Insert loop over spacial coordinate zeta for substrate diffusion,
    % particulates, biomass growth within biofilm, etc
    
    %Call on Biofilm Surface Substrate Concentration from 'Diffusion'
    
    %Calculate Input Parameters
    Sbold=Sb; %Solve Diffusion with previous solution for efficiency
    z=linspace(0,param.Lf,Nz); %[m] Grid of Biofilm Depth
    dz=z(2)-z(1); %[m]
    [Cs,Sb,bflux(i+1),flux(i+1)]=biofilmdiffusion(Sbold,S(i),Nz,dz,param);
    
    %Call on Biofilm Thickness and Vdet/Vg from 'BiofilmThickness_Fn'
    Lf_old=param.Lf;
    [param.Lf,Vdet]=lf(Sb,Lf_old,dt,dz,param);
    
    %Call on tank substrate/biomass concentration from 'tankenvironment'
    [t(i+1),x(i+1),S(i+1)]=tankenvironment(t(i),x(i),S(i),SA,Vdet,dt,Cs,Co,param);
   
    %Call on desired plots from 'outputs'
    outIter=outIter+1;
    if (outIter==outFreq)
        [plots] = outputs(t(1:i+1),x(1:i+1),S(1:i+1),z,bflux(1:i+1),Sb,param.Lf,plots);
        outIter=0;
    end
    
end
%end
toc