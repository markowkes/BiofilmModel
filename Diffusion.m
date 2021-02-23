function [Cs,Sb,bflux]=Diffusion(Sbold,Lf,LL,So,mumax,Xb,Yxs,De,Km,Daq)
%This Function will take initial tank conditions and model the diffusion of
% substrates into the biofilm. The results of this uptake will be used to
% model the manner in which tank conditions reach equilibrium

%Iterations
tic
for iter=1:1000

        c=2:1:Nz-1; %array to run concentrations through
        Sbnew(c)=(Sb(c+1)+Sb(c-1)-(mu(Sb(c),mumax,Km)*Xb*(dz^2))/(Yxs*De))/2; %Concentration of substrate at biofilm depth
    
        %Boundary Conditions After Iteration
        Sbnew(1)=Sbnew(2);
        Sbnew(end)=Sb(end);
        
        %Flux Calculations
              bflux=(Sbnew(end)-Sbnew(end-1))/dz; %Biofilm Flux at boundary
              
        %Flux Matching       
        Sbnew(end)=So-LL*De*(bflux/Daq);
        
        %Non Zero Condition
        Sbnew(Sbnew < 0) = 0;

         Sb=Sbnew;
         
end
Cs=Sb(end); %output Surface Concentration
end
