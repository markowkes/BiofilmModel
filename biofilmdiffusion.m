function [Cs,Sb,bflux]=biofilmdiffusion(Sbold,S,Lf,Nz,dz,LL,mumax,Xb,Yxs,De,Km,Daq)
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
        Sb(c)=(Sbold(c+1)+Sbold(c-1)-(mu(Sbold(c),mumax,Km)*Xb*(dz^2))/(Yxs*De))/2; %Concentration of substrate at biofilm depth
        Sb(c)=lamda*Sb(c)+(1-lamda)*Sbold(c); %Over Relaxation Modification
        
        %Boundary Conditions
        Sb(1)=Sb(2); %Zero Flux at bottom Boundary
        
        if LL<zeroLL
            Sb(end)=S;
        else
            Sb(end)=((Daq/LL)*S+(De/dz)*Sb(end-1))/((De/dz)+(Daq/LL)); %Flux Matching At Top Border
        end
           %Flux Calculations
           bflux=(Sb(end)-Sb(end-1))/dz; %Biofilm Flux at boundary
           Sbold=Sb;
           
        %Non Zero Condition
        Sb(Sb < 0) = 0;
        
          if max(abs(Sb-Sbold))<tol
              break
          end
end
Cs=Sb(end); %output Surface Concentration
end