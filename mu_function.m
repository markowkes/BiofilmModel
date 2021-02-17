function [muSb,muS] = mu_function(mumax,Km,Sb,S)
% This function will allow for m to be calculated from the Monod Growth
% Kinetics Equation
muSb = @(Sb) ((mumax*Sb)/(Km+Sb));
muS = @(S) ((mumax*S)/(Km+S));

end