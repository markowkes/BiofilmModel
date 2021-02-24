function [Cs,Sb,bflux]=Diffusion(Sbold,LL,So,mumax,Xb,Yxs,De,Km,Daq,dz,Nz)
%This Function will take initial tank conditions and model the diffusion of
% substrates into the biofilm. The results of this uptake will be used to
% model the manner in which tank conditions reach equilibrium

Sb=zeros(1,Nz);

%Iterations
for iter=1:1000
muiter=mu(Sbold,mumax,Km); %To compute Sb
        c=2:1:Nz-1; %array to run concentration calculations through
        Sb(c)=(Sbold(c+1)+Sbold(c-1)-(muiter(c)*Xb*(dz^2))/(Yxs*De))/2; %Concentration of substrate at biofilm depth
    
        %Boundary Conditions After Iteration
        Sb(1)=Sb(2);
        Sb(end)=Sbold(end);
        
        %Flux Calculations
        bflux=(Sb(end)-Sb(end-1))/dz; %Biofilm Flux at boundary
              
        %Flux Matching       
        Sb(end)=So-LL*De*(bflux/Daq); %If LL=0 Boundary Concentration is Tank Concentration
        
        %Non Zero Condition
        Sb(Sb < 0) = 0;

         Sbold=Sb;
end
Cs=Sb(end); %output Surface Concentration
end
