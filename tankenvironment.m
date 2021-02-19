function [S,x,t]=tankenvironment(xo,So,V,SA,Qdot,Sin,Vdet,mumax,Km,Yxs,Daq,LL,Cs,Co,Xb,dt,N,i)
%% This function describes the greater tank environment and assumed that it is well mixed
% It calls all the necessary tank geometry, flow parameters, and specific
% parameters describing the biofilm and uses the differential equations
% describing the substrate and biomass concentrations through time to
% produce plots profiling these concetrations over a set period of time

dxdt = @(x,t,S,Cs,mu,Vdet) (mu-(Qdot/V))*x+Vdet*SA*Xb; %Biomass Concentration Change wrt time
dsdt = @(x,t,S,Cs,mu) -((mu*x)/Yxs)+((Qdot*Sin)/V)-((Qdot*S)/V)-(SA*((Daq/LL)*(Co-Cs))); 
% ^^^Substrate Concentration Change wrt time, now also considers flux
% through boundary layer of biofilm

%Preallocation
t = zeros(1,N); %Time
x = zeros(1,N); %Biomass Concentration in bulk liquid
S = zeros(1,N); %Substrate in bulk liquid


%Initial Conditions
t(1)=0;
x(1)=xo;
S(1)=So;


    
t(i+1) = t(i) + dt;
    
xstar = x(i) + dt*dxdt(x(i),t(i),S(i),Cs,mu(S(i),mumax,Km),Vdet);
Sstar = S(i) + dt*dsdt(x(i),t(i),S(i),Cs,mu(S(i),mumax,Km));
    
x(i+1) = x(i) + dt/2*(dxdt(x(i),t(i),S(i),Cs,mu(S(i),mumax,Km),Vdet)+dxdt(xstar,t(i+1),Sstar,Cs,mu(S(i+1),mumax,Km),Vdet));
S(i+1) = S(i) + dt/2*(dsdt(x(i),t(i),S(i),Cs,mu(S(i),mumax,Km)+dsdt(xstar,t(i+1),Sstar,Cs,mu(S(i+1),mumax,Km)))); 
    

end
