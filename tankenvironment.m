function [Snew,xnew,tnew]=tankenvironment(t,x,S,V,SA,Qdot,Sin,Vdet,mumax,Km,Yxs,Daq,LL,Cs,Co,Xb,dt)
%% This function describes the greater tank environment and assumed that it is well mixed
% It calls all the necessary tank geometry, flow parameters, and specific
% parameters describing the biofilm and uses the differential equations
% describing the substrate and biomass concentrations through time to
% produce plots profiling these concetrations over a set period of time

dxdt = @(x,t,S,Cs,Vdet) (mu(S,mumax,Km)-(Qdot/V))*x+Vdet*SA*Xb; %Biomass Concentration Change wrt time
dsdt = @(x,t,S,Cs) -((mu(S,mumax,Km)*x)/Yxs)+((Qdot*Sin)/V)-((Qdot*S)/V)-(SA*((Daq/LL)*(Co-Cs))); 
% ^^^Substrate Concentration Change wrt time, now also considers flux
% through boundary layer of biofilm

    
tnew = t + dt;
    
xstar = x + dt*dxdt(x,t,S,Cs,Vdet);
Sstar = S + dt*dsdt(x,t,S,Cs     );
    
xnew = x + dt/2*(dxdt(x,t,S,Cs,Vdet)+dxdt(xstar,tnew,Sstar,Cs,Vdet));
Snew = S + dt/2*(dsdt(x,t,S,Cs     )+dsdt(xstar,tnew,Sstar,Cs     )); 
    

end
