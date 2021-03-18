function tankenvironment(t,x,S,SA,Vdet,dt,Cs,Co,param)
    ## This function describes the greater tank environment and assumed that it is well mixed
    # It calls all the necessary tank geometry, flow parameters, and specific
    # parameters describing the biofilm and uses the differential equations
    # describing the substrate and biomass concentrations through time to
    # produce plots profiling these concetrations over a set period of time
    Q=param.Q
    V=param.V
    Xb=param.Xb
    Yxs=param.Yxs
    Sin=param.Sin
    Daq=param.Daq
    LL=param.LL
    
    #dxdt = x,S -> (mu(S,param)-(Q/V))*x+Vdet*SA*Xb
    #dsdt = x,S -> -((mu(S,param)*x)/Yxs)+((Q*Sin)/V)-((Q*S)/V)-(SA*((Daq/LL)*(Co-Cs)))
    
    # Packing y
    y=[x, S]
    
    #f = t,y -> [dxdt(y[0],y[1]),
    #            dsdt(y[0],y[1])]

    function f(t,y)
        x=y[1]
        S=y[2]
        rhs=zeros(2)
        rhs[1]=(mu(S,param)-(Q/V))*x+Vdet*SA*Xb
        rhs[2]=-((mu(S,param)*x)/Yxs)+((Q*Sin)/V)-((Q*S)/V)-(SA*((Daq/LL)*(Co-Cs)))
        return rhs
    end
           
    tol=1e-8
    tnew=0
    ynew=zeros(size(y))
    local iter
    for outer iter=1:10
        s1 = f(t     ,y            )
        s2 = f(t+  dt/2,y+  dt/2*s1)
        s3 = f(t+3*dt/4,y+3*dt/4*s2)
        
        tnew = t + dt
        ynew = y + dt/9*(2*s1 + 3*s2 + 4*s3)
        
        s4 = f(tnew,ynew)
        
        error = dt/72*(-5*s1 + 6*s2 + 8*s3 - 9*s4)
        
        # Update timestep
        if all(abs.(error) .< tol/100)
            # dt is getting very small
            dt=dt*2
        elseif any(abs.(error) .> tol)
            # dt is too big
            dt=dt/2
        else
            # Step completed with good dt
            break
        end 
    end 
    
    # Unpacking y
    xnew=ynew[1]
    Snew=ynew[2]
        
    return tnew,xnew,Snew,dt,iter
end