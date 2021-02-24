function [Cs,Sb,bflux]=diffusion(Sbold,LL,S,mumax,Xb,Yxs,De,Km,Daq,Lf,dz)
%% This function models the diffusion of a substrate within the biofilm
%This Function will take tank conditions (So,Xb,LL) and various growth factors (Yxs,De,Km,Daq) model the diffusion of
% substrates into the biofilm over the grid . The results of this uptake will be used to
% model the manner in which tank conditions reach equilibrium
 Nz=50;
 
 Sb=zeros(1,Nz);

%Iterations
for iter=1:1000
%muiter=mu(Sb,mumax,Km); %To compute Sb
        c=2:1:Nz-1; %array to run concentration calculations through
        Sb(c)=(Sbold(c+1)+Sbold(c-1)-(mu(Sbold(c),mumax,Km)*Xb*(dz^2))/(Yxs*De))/2; %Concentration of substrate at biofilm depth
    
        %Boundary Conditions After Iteration
        Sb(1)=Sb(2);
        Sb(end)=S;

        %Flux Calculations
        bflux=(Sb(end)-Sb(end-1)); %Biofilm Flux at boundary
              
        %Flux Matching       
        %Sb(end)=S-LL*De*(bflux/(dz*Daq)); %If LL=0 Boundary Concentration is Tank Concentration
        
         %Non Zero Condition
        Sb(Sb < 0) = 0;
        
         Sbold=Sb;
end
Cs=Sb(end); %output Surface Concentration
end