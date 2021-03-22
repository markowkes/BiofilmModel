 %MAINDRIVER
clear; clc
tic
%% Rewriting inputs section to run through test cases
%Call on specific test case parameters
num=7 ; %number of case, A corresponds to 1, B corresponds to 2....
param=cases(num); %structure variables are stored in

%Create initial biofilm grid
Nz=50; %Linear GridPoints in Biofilm
z=linspace(0,param.Lfo,Nz); %[m] Grid of Biofilm Depth
dz=z(2)-z(1); %[m]

%Initial Boundary Conditions (in Biofilm)
Sb=zeros(1,Nz);
Sb(end)=param.So; %initially assume boundary concentration = So

%Time Constraints
tFin=30; %[days]
dt=1e-2; %Interval
N=tFin/dt; %Number of steps
outFreq=2000; %Number of steps between plot updates.

%Preallocation
t       =zeros(1,N); %Time
x       =zeros(1,N); %Biomass Concentration in bulk liquid
S       =zeros(1,N); %Substrate in bulk liquid
bflux   =zeros(1,N); %Boundary Layer Flux of Biofilm Preallocate
flux    =zeros(1,N); %Right hand side of power point equation to ensure matching flux
Lf      =zeros(1,N); %Right hand side of power point equation to ensure matching flux

%Initial Conditions
t(1)=0;
x(1)=param.xo;
S(1)=param.So;
Lf(1)=param.Lfo;

%Initialize plots 
outIter=outFreq-1;
plots=0; titles=0;

%% Time Loop
i=1;
while t(i)<tFin-dt
    
    %Update biofilm grid as biofilm grows
    z=linspace(0,Lf(i),Nz); %[m] Grid of Biofilm Depth
    dz=z(2)-z(1); %[m]
    
    %Call on "biofilmdiffusion"
    [Sb,bflux(i+1)]=biofilmdiffusion_fd(Sb,S(i),Nz,dz,t(i),param);
    
    %Call on "lf"
    [Lf(i+1),Vdet]=lf(Sb,Lf(i),dt,dz,param);
    
    %Call on "tankenvironment"
    [t(i+1),x(i+1),S(i+1),dt]=tankenvironment(t(i),x(i),S(i),Vdet,dt,bflux(i+1),param);
    
    %Call on desired plots from 'outputs'
    outIter=outIter+1;
    if (outIter>=outFreq)
        [plots,titles] = outputs(t(1:i+1),x(1:i+1),S(1:i+1),z,bflux(1:i+1),Lf(1:i+1),Sb,param,plots,titles);
        outIter=0;
    end
    
    % Update iterator
    i=i+1;
end

% Make final figures
[plots,titles] = outputs(t(1:i),x(1:i),S(1:i),z,bflux(1:i),Lf(1:i),Sb,param,plots,titles);

toc