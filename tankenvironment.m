function [s4,tnew,xnew,Snew,dt]=tankenvironment(t,x,S,Vdet,dt,bflux,param)
%% This function describes the greater tank environment and assumed that it is well mixed
% It calls all the necessary tank geometry, flow parameters, and specific
% parameters describing the biofilm and uses the differential equations
% describing the substrate and biomass concentrations through time to
% produce plots profiling these concetrations over a set period of time

Q  =param.Q;
V  =param.V;
Xb =param.Xb;
Yxs=param.Yxs;
Sin=param.Sin;
A  =param.A;
tol=param.ttol;
dtmax=param.dtmax;

dxdt = @(x,t,S,Vdet) (mu(1,S(:,1),param)-(Q/V))*x+(Vdet*A*Xb)/V; %Biomass Concentration Change wrt time
dsdt = @(x,t,S,k) -((mu(1,S(:,1),param)*x)./Yxs)+((Q.*Sin(:,1))./V)-((Q.*S(:,1))./V)-((A.*bflux)./V); % ^^^Substrate Concentration Change wrt time

% Packing y
y=[x; S(:,1)];

f =@(t,y) [dxdt(y(1),t,y,Vdet)
           dsdt(y(1),t,y(2:end),[1:2])];

while true
%     s1 = f(t     ,y(1:2,:)            );
%     s2 = f(t+  dt/2,y(1:2)+  dt/2*s1);
%     s3 = f(t+3*dt/4,y(1:2)+3*dt/4*s2);
%     
%     tnew = t + dt;
%     ynew = y + dt/9*(2*s1 + 3*s2 + 4*s3);
%     
%     s4 = f(tnew,ynew);
%     
%     error = dt/72*(-5*s1 + 6*s2 + 8*s3 - 9*s4);
    
    s1 = f(t     ,y            );
    s2 = f(t+  dt/2,y+  dt/2*s1);
    s3 = f(t+3*dt/4,y+3*dt/4*s2);
    
    tnew = t + dt;
    ynew = y + dt/9*(2*s1 + 3*s2 + 4*s3);
    
    s4 = f(tnew,ynew);
    
    error = dt/72*(-5*s1 + 6*s2 + 8*s3 - 9*s4);
    
%     dt=1/(0.05*norm(s4));
%     if norm(s4)==0
%         dt=1;
%     else
%         dt=dt;
%     end
    
  
    % Update timestep
    if abs(error) < tol/100 
        %dt is getting very small
        dt=dt*2;
        %Check if dt max is exceded
        if dt>dtmax 
            %Set dt to dtmax
            dt=dtmax;
            break
        end
    elseif abs(error) > tol
        %dt is too big
        dt=dt/2;
    else
        %Step completed with good dt
        break
    end
        
end

% Unpacking y
xnew=ynew(1);
Snew=ynew(2:end);

end
