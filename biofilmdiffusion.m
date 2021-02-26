function [Cs,Sb,bflux]=biofilmdiffusion(Sbold,LL,S,mumax,Xb,Yxs,De,Km,Daq,Lf,dz)
%% This function models the diffusion of a substrate within the biofilm
%This Function will take tank conditions (So,Xb,LL) and various growth factors (Yxs,De,Km,Daq) model the diffusion of
% substrates into the biofilm over the grid . The results of this uptake will be used to
% model the manner in which tank conditions reach equilibrium
Nz=50;
Sb=zeros(1,Nz);

 
lamda=.5; %Factor For Over Relaxation Method
tol=1e-2;

%Iterations
for iter=1:10000
        c=2:1:Nz-1; %array to run concentration calculations through
        Sb(c)=(Sbold(c+1)+Sbold(c-1)-(mu(Sbold(c),mumax,Km)*Xb*(dz^2))/(Yxs*De))/2; %Concentration of substrate at biofilm depth
        Sb(c)=lamda*Sb(c)+(1-lamda)*Sbold(c); %Over Relaxation Method
        
        %Boundary Conditions After Iteration
        Sb(1)=Sb(2); %Zero Flux at bottom Boundary
%        Sb(end)=S;
%        Sb(end)=Sbold(end);
        Sb(end)=((Daq/LL)*S+(De/dz)*Sb(end-1))/((De/dz)+(Daq/LL)); %Flux Matching At Top Border

%          %Flux Matching     
%          Sb(end)=S-LL*De*(bflux/(Daq)); %If LL=0 Boundary Concentration is Tank Concentration

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