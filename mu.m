function [mu] = mu(S,mumax,Km)
% This function will allow for m to be calculated from the Monod Growth
% Kinetics Equation
mu = ((mumax*S)./(Km+S));

end