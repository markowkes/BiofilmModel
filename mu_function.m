function [mu] = mu_function(mumax,Km,S)
% This function will allow for m to be calculated from the Monod Growth
% Kinetics Equation
mu = @(S) ((mumax*S)/(Km+S)); 
end