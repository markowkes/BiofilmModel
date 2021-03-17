function diffusion(Sbold,S,Nz,dz,param)
    
    tol=1e-12
    
    # Get variables out out of param
    Xb=param.Xb
    Yxs=param.Yxs
    De=param.De
    LL=param.LL
    Daq=param.Daq
    
    # RHS of eqn for interior points
    g = S -> mu(S,param)*Xb/(Yxs*De)
    
    # Initial guess
    Sb=Sbold
    
    # Preallocate matrices
    A=zeros(Nz,Nz)
    for i=2:Nz-1
        A[i,i-1]=-1 # Lower diagonal
        A[i,i+1]=-1 # Upper diagonal
    end
    B=zeros(Nz)

    # BC at top of biofilm
    A[Nz,Nz  ]= De*LL+Daq*dz
    A[Nz,Nz-1]=-De*LL
    B[Nz]     = Daq*dz*S
    
    # BC at bottom of biofilm 
    A[1,1]= 1
    A[1,2]=-1
    B[1]  = 0
    
    for iter=1:100
        
        # Interior Points
        delta=1e-3
        Sb_p=Sb[1:Nz-1].+delta
        Sb_m=Sb[1:Nz-1].-delta
        dgds=(g(Sb_p)-g(Sb_m))/(Sb_p-Sb_m)
        for i=2:Nz-1
            A[i,i]=2+dz^2*dgds[i-1]
            B[i]=dz^2*(Sb[i]*dgds[i]-g(Sb[i]))
        end

        # Solve for new concentration
        Sb = A\B
        
        # Non Zero Condition
        Sb=max.(0,Sb)

        # Check if converged
        if maximum(abs.(Sb-Sbold)) < tol
            break
        end
    
        # Transfer solution for next iteration
        Sbold=Sb
    end

    Cs=Sb[Nz] # output Surface Concentration

    # Flux Calculations
    bflux=De *(Sb[Nz]-Sb[Nz-1])/dz # Biofilm Flux
    flux =Daq*(S     -Sb[Nz  ])/LL # Boundary Layer Flux
            
    return Cs,Sb,bflux,flux
end