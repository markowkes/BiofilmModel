function [t,x,S,Lf]=MAINDRIVER(param)

Nz      =param.Nz;
Nx      =param.Nx;
Ns      =param.Ns;
So      =param.So;
phio    =param.phio;
rho     =param.rho;
tFin    =param.tFin;
outFreq =param.outFreq;

%% Preallocation
%Compute number of time-steps to solve for
N=round(tFin/param.dtmax); 

%Corresponding arrays
t       =zeros(1,N); %Time
x       =zeros(Nx,N); %Biomass Concentration in bulk liquid
S       =zeros(Ns,N); %Substrate in bulk liquid
bflux   =zeros(Ns,N); %Boundary Layer Flux of Biofilm Preallocate
Lf      =zeros(1,N); %Right hand side of power point equation to ensure matching flux
dt      =zeros(1,N); %size of each time step

%% Initial Conditions
%Biofilm
Xb=zeros(Nx,Nz);
phi=zeros(Nx,Nz);
for i=1:Nx
    phi(i,:)=phio(i);
end
for i=1:Nz % Place into own function
    for j=1:Nx
        Xb(j,i)=rho(j).*phi(j,i);
    end
end
Sb=zeros(Ns,Nz);
Sb(:,end)=So; %initially assume boundary concentration = So 

Lf(1)=param.Lfo;

%Tank
t(1)=0;
x(:,1)=param.xo;
S(:,1)=So;

%Time
dt(1)=param.dtmax;

%% Initialize plots 
outIter=outFreq-1;
plots=0; titles=0;

%% Time Loop
i=1;
while t(i)<tFin-dt(i)
    
    % Check if arrays are filling up
    if length(bflux)==i 
        % Compute an estimate for reamaining time steps
        Nrem    =round((tFin-t(i))/dt(i));
        
        % Append time dependant arrays with estimate
        t       =[t     zeros(1,Nrem)]; 
        x       =[x     zeros(Nx,Nrem)]; 
        S       =[S     zeros(Ns,Nrem)];
        bflux   =[bflux zeros(Nx,Nrem)]; 
        Lf      =[Lf    zeros(1,Nrem)]; 
        dt      =[dt    zeros(1,Nrem)];        
    end
    
    %Update biofilm grid as biofilm grows
    z=linspace(0,Lf(i),Nz); %[m] Grid of Biofilm Depth
    dz=z(2)-z(1); %[m]
    
    %Call on "biofilmdiffusion"
    [Sb,bflux(:,i+1)]=biofilmdiffusion_fd(Sb,S(:,i),Xb,dz,t(i),param);
    
    %Call on "particulates"
    [Lf(i+1),Vdet,phi,Xb]=particulates(Sb,phi,Lf(i),dt(i),dz,param);
    
    %Call on "tankenvironment"
    [t(i+1),x(:,i+1),S(:,i+1),dt(i+1)]=tankenvironment(t(i),x(:,i),S(:,i),Vdet,Xb,dt(i),bflux(:,i+1),param);
    
    %Call on desired plots from 'outputs'
    outIter=outIter+1;
    if (outIter>=outFreq)
        [plots,titles] = outputs(t(1:i+1),x(:,1:i+1),S(:,1:i+1),z,bflux(:,1:i+1),Lf(1:i+1),Sb(:,:),param,plots,titles);
        outIter=0;
    end
    
    % Update iterator
    i=i+1;
end

% Make final figures
[~,~] = outputs(t(1:i),x(:,1:i),S(:,1:i),z,bflux(:,1:i),Lf(1:i),Sb,param,plots,titles);

% Remove extra zeros if they exist
t=t(1:i);
x=x(:,1:i);
S=S(:,1:i);
Lf=Lf(1:i);

end
