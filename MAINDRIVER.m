function [t,x,S,Lf]=MAINDRIVER(param)

Xbo     =param.Xb;
Nz      =param.Nz;
Nx      =param.Nx;
Ns      =param.Ns;
So      =param.So;
phio    =param.phio;
tFin    =param.tFin;
outFreq =param.outFreq;

%% Preallocation
%Compute number of time-steps to solve for
N=round(tFin/param.dtmax); 

%Corresponding arrays
t       =zeros(1,N); %Time
x       =zeros(Nx,N); %Biomass Concentration in bulk liquid
S       =zeros(2,N); %Substrate in bulk liquid
bflux   =zeros(2,N); %Boundary Layer Flux of Biofilm Preallocate
Lf      =zeros(1,N); %Right hand side of power point equation to ensure matching flux
dt      =zeros(1,N); %size of each time step

%% Initial Conditions
%Biofilm
Xb=zeros(2,Nz);
Xb(1,:)=Xbo(1);
Xb(2,:)=Xbo(2);
Sb=zeros(2,Nz);
phi=zeros(2,Nz);
phi(1,:)=phio(1);
phi(2,:)=phio(2);
Sb(:,end)=So; %initially assume boundary concentration = So 

Lf(1)=param.Lfo;

%Tank
t(1)=0;
x(1)=param.xo;
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
    [Lf(i+1),Vdet,Xb]=particulates(Sb,phi,Lf(i),dt(i),dz,param);
    
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
