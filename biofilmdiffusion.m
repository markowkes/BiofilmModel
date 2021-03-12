function [Cs,Sb,bflux,flux]=biofilmdiffusion(Sbold,S,Nz,dz,param)
%% This function models the diffusion of a substrate within the biofilm
%This Function will take tank conditions (So,Xb,LL) and various growth factors (Yxs,De,Km,Daq) model the diffusion of
% substrates into the biofilm over the grid . The results of this uptake will be used to
% model the manner in which tank conditions reach equilibrium

Sb=Sbold; %preallocate array
lamda=.2; %Factor For Over Relaxation Method

zeroLL=1e-10; %[m] condition to consider zero thickness boundary layer
tol=1e-8; %tolerance for conversion

% Get variables out out of param
Xb=param.Xb;
Yxs=param.Yxs;
De=param.De;
LL=param.LL;
Daq=param.Daq;

%Iterations
for iter=1:100000
    c=2:1:Nz-1; %array to run concentration calculations through
    Sb(c)=(Sbold(c+1)+Sbold(c-1)-(mu(Sbold(c),param)*Xb*(dz^2))/(Yxs*De))/2; %Concentration of substrate at biofilm depth
    Sb(c)=lamda*Sb(c)+(1-lamda)*Sbold(c); %Over Relaxation Modification
        
    %Boundary Conditions
    Sb(1)=Sb(2); %Zero Flux at bottom Boundary
    % Top of biofilm
    if LL>zeroLL
        % Flux Matching
        Sb(end)=((Daq/LL)*S+(De/dz)*Sb(end-1))/((De/dz)+(Daq/LL)); 
    else
        % Set concentration at top to S
        Sb(end)=S;
    end
                 
    %Non Zero Condition
    Sb(Sb < 0) = 0;
    
    % Check if converged
    if max(abs(Sb-Sbold))<tol
        break
    end

    % Transfer solution for next iteration
    Sbold=Sb;
end
Cs=Sb(end); %output Surface Concentration

%Flux Calculations
bflux=De *(Sb(end)-Sb(end-1))/dz; %"Biofilm Flux" at boundary LHS of provided flux matching equation
flux =Daq*(S      -Sb(end  ))/LL; %"Boundary Layer Flux" RHS of provided flux matching equation
    
end