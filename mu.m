function [mu] = mu(S,param)
%% This function will allow for mu to be calculated 
% This function calls various growth models to calculate an
% array of growth rates mu for any array 'S'
mumax=param.mumax;
Km=param.Km;

mu = ((mumax*S)./(Km+S));                      % Monod Growth Rate Equation
% mu2 = ((mumax*Sa)./(Kma+Sa))*(Sb./(Kmb+Sb));    % Double Monod Growth Rate Equation
% mu3 = ((mumax*Sa)./(Kma+Sa))*(1./(1+(Sb/Kmb))); % Inhibition Growth Rate Equation
% mu4 = 0;                                        % None                 
end