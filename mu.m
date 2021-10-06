function [mu] = mu(j,S,param)
%% This function will allow for mu to be calculated 
% This function calls various growth models to calculate an
% array of growth rates mu for any array 'S'
model=param.model;
mumax=param.mumax;
Km=param.Km;

switch model(j)
    case 1 %'Linear Growth Rate Equation'
        mu = ((mumax*S(j))./(Km(1,1,j)));
    case 2 %'Monod Growth Rate Equation'
        mu = ((mumax*S(j))./(Km(1,1,j)+S(j)));   
    case 3 %'Double Monod Growth Rate Equation'
        mu=1;
        for i=1:length(S)
            mu = mu*(mumax*S(i))./(Km(2,2,i)+S(i))
        end          
%         Sa=S(j,:);
%         Sb=S(Km(1,2,j),:);
% 
%         Kma=Km(2,1,j);
%         Kmb=Km(2,2,j);
%         
%         mu = ((mumax*Sa)./(Kma+Sa)).*(Sb./(Kmb+Sb));
    case 4 %'Inhibition Growth Rate Equation'
        
        Sa=S(1,:);
        Sb=S(2,:);

        Kma=Km(1,1,1);
        Kmb=Km(1,1,2);
        
        mu = ((mumax*Sa)./(Kma+Sa))*(1./(1+(Sb/Kmb)));
    case 5 %'None'
        mu = 0;
    otherwise
        warning('unrecognized Growth Rate Model specified in "cases"')
end
end