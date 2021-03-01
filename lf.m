function [Lf,Vdet]=lf(Sb,Lf_old,dt,dz,param)
%This Function takes the substrate concentration at a given instant in time
%(and the corresponding equation for mu) as well as the old Biofilm
%thickness and computes a new Biofilm Thickness (Lf)

%Corresponding outputs for the growth and detachment rates are computed for
%a given instant of time.

% Compute mean mu - growthrate
muBar=mean(mu(Sb(:),param.mumax,param.Km));

%Biofilm Thickness
Lf=Lf_old+dt*(muBar*Lf_old-param.Kdet*Lf_old^2); %New Biofilm thickness at instant

%Detachment
Vdet=param.Kdet*Lf^2; %New %Velocity of mass leaving biofilm into bulk liquid

%Growth
Vg=0; %initial condition for loop
for i=1:length(Sb)-1
    Vg=Vg+dz*((mu(Sb(i),param.mumax,param.Km)+mu(Sb(i+1),param.mumax,param.Km))/2); %trapezoidal integration method for
                                                            %new growth velocity of biofilm
end
end
