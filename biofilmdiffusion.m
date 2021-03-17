function [Cs,Sb,bflux,flux]=biofilmdiffusion(Sbold,S,Nz,dz,t,param)
%% This function models the diffusion of a substrate within the biofilm
%This Function will take tank conditions (So,Xb,LL) and various growth factors (Yxs,De,Km,Daq) model the diffusion of
% substrates into the biofilm over the grid . The results of this uptake will be used to
% model the manner in which tank conditions reach equilibrium

Sb=Sbold; %preallocate array
lamda=.2; %Factor For Over Relaxation Method

tol=1e-8; %tolerance for conversion

% Get variables out out of param
Xb=param.Xb;
Yxs=param.Yxs;
De=param.De;
LL=param.LL;
Daq=param.Daq;

%Iterations
iter=10000;

%Solve
for i=1:iter
    c=2:1:Nz-1; %array to run concentration calculations through
    Sb(c)=(Sbold(c+1)+Sbold(c-1)-(mu(Sbold(c),param)*Xb*(dz^2))/(Yxs*De))/2; %Concentration of substrate at biofilm depth
    Sb(c)=lamda*Sb(c)+(1-lamda)*Sbold(c); %Over Relaxation Modification
        
    % Boundary Conditions
    Sb(1)=Sb(2); %Zero Flux at bottom Boundary
    % Top of biofilm
    
    % Flux Matching
    Sb(end)=(dz*Daq*S+LL*De*Sb(end-1))/(LL*De+dz*Daq); 
                 
    % Non Zero Condition
    Sb(Sb < 0) = 0;
    
    % Check if converged
    if max(abs(Sb-Sbold))<tol
        break
    else
        if i==iter
             fprintf('Diffusion Unable to Converge at time %3.8f\n',t)
        end
    end

    % Transfer solution for next iteration
    Sbold=Sb;
end
Cs=Sb(end); % Output Surface Concentration

% Flux Calculations
bflux=De *(Sb(end)-Sb(end-1))/dz; %"Biofilm Flux" at boundary LHS of provided flux matching equation
flux =Daq*(S      -Sb(end  ))/LL; %"Boundary Layer Flux" RHS of provided flux matching equation
    
end