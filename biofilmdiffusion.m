function [Cs,Sb,bflux,flux]=biofilmdiffusion(Sbold,S,Nz,dz,param)
%% This function models the diffusion of a substrate within the biofilm
%This Function will take tank conditions (So,Xb,LL) and various growth factors (Yxs,De,Km,Daq) model the diffusion of
% substrates into the biofilm over the grid . The results of this uptake will be used to
% model the manner in which tank conditions reach equilibrium

Sb=zeros(1,Nz); %preallocate array
lamda=.5; %Factor For Over Relaxation Method

zeroLL=1e-10; %[m] condition to consider zero thickness boundary layer
tol=1e-2; %tolerance for conversion

%Iterations
for iter=1:10000
        c=2:1:Nz-1; %array to run concentration calculations through
        Sb(c)=(Sbold(c+1)+Sbold(c-1)-(mu(Sbold(c),param)*param.Xb*(dz^2))/(param.Yxs*param.De))/2; %Concentration of substrate at biofilm depth
        Sb(c)=lamda*Sb(c)+(1-lamda)*Sbold(c); %Over Relaxation Modification
        
        %Boundary Conditions
        Sb(1)=Sb(2); %Zero Flux at bottom Boundary
        
%         if LL<zeroLL
%             Sb(end)=S;
%         else
            Sb(end)=((param.Daq/param.LL)*S+(param.De/dz)*Sb(end-1))/((param.De/dz)+(param.Daq/param.LL)); %Flux Matching At Top Border
%         end
           %Flux Calculations
           bflux=(Sb(end)-Sb(end-1))*(1/dz); %"Biofilm Flux" at boundary LHS of provided flux matching equation
           flux=(S-Sb(end))*(param.Daq/param.LL); %"Boundary Layer Flux" RHS of provided flux matching equation
           Sbold=Sb;
           
        %Non Zero Condition
        Sb(Sb < 0) = 0;
        
          if max(abs(Sb-Sbold))<tol
              break
          end
end
Cs=Sb(end); %output Surface Concentration
end