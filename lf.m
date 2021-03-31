function [Lf,Vdet]=lf(Sb,Lf_old,dt,dz,param)
%This function takes the substrate concentration at a given instant in
%time and the old biofilm thickness to computes the growth velocity as well
%as the detachement velocity at a given instant in time. These results are
%used to compute a new biofilm thickness Lf

%Compute mean mu - growthrate
%muBar=mean(mu(Sb(:),param));

Kdet=param.Kdet;

%% Growth

% %Trapzoidal Integration Method
Vg=0; %initial growth velocity
for i=1:length(Sb)-1
    Vg=Vg+dz*((mu(Sb(i),param)+mu(Sb(i+1),param))/2);
end

% %Midpoint Integration Method
% Vg=0; %initial growth velocity
% i=1;
% while i<=length(Sb)-2
%     Vg=Vg+dz*(mu(Sb(i),param)+mu(Sb(i+2),param));
%     i=i+2;
% end

% %Simpsons Integration
% Vg=0;
% if rem(length(Sb),2)~=1
%     fprintf("grid size needs to be odd!\n")
% end
% Vg=(dz/3)*(mu(Sb(1),param) + mu(Sb(end),param) ...
%    + 4*sum(mu(Sb(2:2:end-1),param)) ...
%    + 2*sum(mu(Sb(3:2:end-2),param)));


%% Detachment
Vdet=Kdet*Lf_old^2; %New Velocity of mass leaving biofilm into bulk liquid

%% Biofilm Thickness
Lf=Lf_old+dt*(Vg-Vdet); %New Biofilm thickness at instant

end
