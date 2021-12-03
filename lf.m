function [Lf,Vdet,phi]=lf(Sb,phi_old,Lf_old,dt,dz,param)
%This function takes the substrate concentration at a given instant in
%time and the old biofilm thickness to computes the growth velocity as well
%as the detachement velocity at a given instant in time. These results are
%used to compute a new biofilm thickness Lf

Kdet=param.Kdet;

%% Growth

Ns       =size(Sb,1);
Nb       =param.Nb;
phi_tot  =param.phi_tot;
dphidt   =zeros(Nb,length(Sb));
phi      =phi_old;

% Compute and store growth velocities
Vg = zeros(1,length(Sb));
for j=1:Nb
    for i=2:length(Sb)-1 % Input boundary condition @ i=1 v=0
        Vg(i)=Vg(i-1)+1/phi_tot*phi(i)*dz*(mu(j,Sb(:,i),param)+mu(j,Sb(:,i-1),param))/2;
    end
end

for j=1:Nb
    for i=2:length(Sb)-1 % Volume fraction boundary condition (Nuemann?)
        dphidt(j,i) = mu(j,Sb(:,i),param)*phi(j,i)-(Vg(i+1)*phi(j,i+1)-Vg(i-1)*phi(j,i-1))/(2*dz);
    end
    i=1; dphidt(j,i) = mu(j,Sb(:,i),param)*phi(j,i);
    i=length(Sb); dphidt(j,i) = mu(j,Sb(:,length(Sb)),param)*phi(j,i)-(Vg(i)*phi(j,i)-Vg(i-1)*phi(j,i-1))/dz;
    % add Nz
end

phi = phi_old+dt*dphidt;

%% Old Integration Methods
% %Trapzoidal Integration Method
% for j=1:Ns
%     for i=1:length(Sb)-1
%         Vg=Vg+dz*((mu(1,Sb(:,i),param)+mu(1,Sb(:,i+1),param))/2);
%     end
% end

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

%% Surface Velocity
Vs = Vg(end)-Vdet;

%% Biofilm Thickness
Lf=Lf_old+dt*(Vs); %New Biofilm thickness at instant??????

end