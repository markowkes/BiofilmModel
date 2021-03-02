function [mu] = mu(S,param)
%% This function will allow for mu to be calculated 
% This function calls on the Monod Growth Kinetics Equation to calculate an
% array of growth rates for any array 'S'
mu = ((param.mumax*S)./(param.Km+S));

end