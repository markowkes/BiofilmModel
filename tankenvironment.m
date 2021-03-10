function [tnew,xnew,Snew]=tankenvironment(t,x,S,SA,Vdet,dt,Cs,Co,param)
%% This function describes the greater tank environment and assumed that it is well mixed
% It calls all the necessary tank geometry, flow parameters, and specific
% parameters describing the biofilm and uses the differential equations
% describing the substrate and biomass concentrations through time to
% produce plots profiling these concetrations over a set period of time

dxdt = @(x,t,S,Cs,Vdet) (mu(S,param)-(param.Q/param.V))*x+Vdet*SA*param.Xb; %Biomass Concentration Change wrt time
dsdt = @(x,t,S,Cs) -((mu(S,param)*x)/param.Yxs)+((param.Q*param.Sin)/param.V)-((param.Q*S)/param.V)-(SA*((param.Daq/param.LL)*(Co-Cs))); 
% ^^^Substrate Concentration Change wrt time, now also considers flux
% through boundary layer of biofilm

    
tnew = t + dt;
    
xstar = x + dt*dxdt(x,t,S,Cs,Vdet);
Sstar = S + dt*dsdt(x,t,S,Cs     );
    
xnew = x + dt/2*(dxdt(x,t,S,Cs,Vdet)+dxdt(xstar,tnew,Sstar,Cs,Vdet));
Snew = S + dt/2*(dsdt(x,t,S,Cs     )+dsdt(xstar,tnew,Sstar,Cs     )); 
    

% dxdt=(mu(S,param)-(param.Q/param.V))*x+Vdet*SA*param.Xb; %Biomass Concentration Change wrt time
% dsdt=-((mu(S,param)*x)/param.Yxs)+((param.Q*param.Sin)/param.V)-((param.Q*S)/param.V)-(SA*((param.Daq/param.LL)*(Co-Cs))); 
% % ^^^Substrate Concentration Change wrt time
% 
% tspan=[0,t];
% [tnew,xnew]=ode15s(dxdt,tspan,x);
% [tnew,Snew]=ode15s(dsdt,tspan,S);
end
