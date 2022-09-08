function [tnew,xnew,Snew,dt]=tankenvironment(t,x,S,Vdet,Xb,dt,bflux,param)
%% This function describes the greater tank environment and assumed that it is well mixed
% It calls all the necessary tank geometry, flow parameters, and specific
% parameters describing the biofilm and uses the differential equations
% describing the substrate and biomass concentrations through time to
% produce plots profiling these concetrations over a set period of time

Q  =param.Q;
V  =param.V;

Yxs=param.Yxs;
Sin=param.Sin;
A  =param.A;
tol=param.ttol;
dtmax=param.dtmax;

Nx = 1;


while true
    % RHS 1
    k1_x = dxdt(t     ,Nx,x,S,Vdet,Q,V,A,Xb,Yxs,param);
    k1_S = dSdt(t     ,Nx,x,S,Q,V,A,Xb,Yxs,Sin,bflux,param);
    
    % RHS 2
    % variables at dt/2
    t2= t+dt/2;
    x2= x+dt/2*k1_x;
    S2= S+dt/2*k1_S;
    k2_x = dxdt(t2,Nx,x2,S2,Vdet,Q,V,A,Xb,Yxs,param);
    k2_S = dSdt(t2,Nx,x2,S2,Q,V,A,Xb,Yxs,Sin,bflux,param);
    
    % RHS 3
    % variables at 3*dt/4
    t3= t+3*dt/4;
    x3= x+3*dt/4*k2_x;
    S3= S+3*dt/4*k2_S;
    k3_x = dxdt(t3,Nx,x3,S3,Vdet,Q,V,A,Xb,Yxs,param);
    k3_S = dSdt(t3,Nx,x3,S3,Q,V,A,Xb,Yxs,Sin,bflux,param);
    
    % Update with computed RHS's
    
    tnew = t + dt;
    xnew = x + dt/9*(2*k1_x + 3*k2_x + 4*k3_x);
    Snew = S + dt/9*(2*k1_S + 3*k2_S + 4*k3_S);
    
    % RHS 4 with updated variables
    k4_x = dxdt(tnew,Nx,xnew,Snew,Vdet,Q,V,A,Xb,Yxs,param);
    k4_S = dSdt(tnew,Nx,xnew,Snew,Q,V,A,Xb,Yxs,Sin,bflux,param);

    error_x = dt/72*(-5*k1_x + 6*k2_x + 8*k3_x - 9*k4_x);
    error_S = dt/72*(-5*k1_S + 6*k2_S + 8*k3_S - 9*k4_S);

    error=max(error_x,error_S);
    
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

end

function dxdt = dxdt(t,Nx,x,S,Vdet,Q,V,A,Xb,Yxs,param) 
    dxdt = zeros(Nx);
    for j=1:Nx
        dxdt(j) = (param.mu{j}(S(:,1),param)-(Q/V))*x(j)+(Vdet(j)*A*Xb(j,end))/V;
    end
end

function dSdt = dSdt(t,Nx,x,S,Q,V,A,Xb,Yxs,Sin,bflux,param) % ^^^Substrate Concentration Change wrt time
    dSdt = ((Q.*Sin(:,1))./V)-((Q.*S(:,1))./V)-((A.*bflux)./V);
    for j=1:Nx
            dSdt = dSdt-((param.mu{j}(S(:,1),param)*x(j))./Yxs(j));
    end
end
%dxbdt= @
