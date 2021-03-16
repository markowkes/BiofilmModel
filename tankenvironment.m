function [tnew,xnew,Snew,dt]=tankenvironment(t,x,S,SA,Vdet,dt,Cs,Co,param)
%% This function describes the greater tank environment and assumed that it is well mixed
% It calls all the necessary tank geometry, flow parameters, and specific
% parameters describing the biofilm and uses the differential equations
% describing the substrate and biomass concentrations through time to
% produce plots profiling these concetrations over a set period of time

Q=param.Q;
V=param.V;
Xb=param.Xb;
Yxs=param.Yxs;
Sin=param.Sin;
Daq=param.Daq;
LL=param.LL;

% dxdt = @(x,t,S,Cs,Vdet) (mu(S,param)-(Q/V))*x+Vdet*SA*Xb; %Biomass Concentration Change wrt time
% dsdt = @(x,t,S,Cs) -((mu(S,param)*x)/Yxs)+((Q*Sin)/V)-((Q*S)/V)-(SA*((Daq/LL)*(Co-Cs))); % ^^^Substrate Concentration Change wrt time
% 
% % Packing y
% y=[x; S];
% 
% f =@(t,y) [dxdt(y(1),t,y(2),Cs,Vdet)
%            dsdt(y(1),t,y(2),Cs)];
%        
% tol=1e-8;
% while true
%     s1 = f(t     ,y            );
%     s2 = f(t+  dt/2,y+  dt/2*s1);
%     s3 = f(t+3*dt/4,y+3*dt/4*s2);
%     
%     tnew = t + dt;
%     ynew = y + dt/9*(2*s1 + 3*s2 + 4*s3);
%     
%     s4 = f(tnew,ynew);
%     
%     error = dt/72*(-5*s1 + 6*s2 + 8*s3 - 9*s4);
%     
%     % Update timestep
%     if abs(error) < tol/100 
%         % dt is getting very small
%         dt=dt*2;
%     elseif abs(error) > tol
%         % dt is too big
%         dt=dt/2;
%     else
%         % Step completed with good dt
%         break
%     end
% end
% 
% % Unpacking y
% xnew=ynew(1);
% Snew=ynew(2);

dxdt = @(x,t,S,Cs,Vdet) (mu(S,param)-(Q/V))*x+Vdet*SA*Xb; %Biomass Concentration Change wrt time
dsdt = @(x,t,S,Cs) -((mu(S,param)*x)/Yxs)+((Q*Sin)/V)-((Q*S)/V)-(SA*((Daq/LL)*(Co-Cs))); %^^^Substrate Concentration Change wrt time

tnew = (t+dt)-0.0001*(abs(gradient(x))/max(norm(x)));
if tnew < (t+dt)/100
    %time step too small
    tnew=tnew*2;
end
%tnew = (t+dt)-0.0001;

xstar = x + dt*dxdt(x,t,S,Cs,Vdet);
Sstar = S + dt*dsdt(x,t,S,Cs     );
    
xnew = x + dt/2*(dxdt(x,t,S,Cs,Vdet)+dxdt(xstar,tnew,Sstar,Cs,Vdet));
Snew = S + dt/2*(dsdt(x,t,S,Cs     )+dsdt(xstar,tnew,Sstar,Cs     ));     


end
