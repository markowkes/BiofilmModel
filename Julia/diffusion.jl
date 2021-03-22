using LinearAlgebra

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
    
    # Diagonals
    Dm=zeros(Nz)   # Main
    Dl=zeros(Nz-1) # Lower  - starts on second row!
    Du=zeros(Nz-1) # Upper

    # RHS
    B=zeros(Nz)

    # Lower/Upper diagonals for interior points
    for i=2:Nz-1
        Dl[i-1]=-1.0  # Dl starts on second row!
        Du[i  ]=-1.0
    end
    
    # Top BC
    Dm[Nz  ]=De*LL+Daq*dz
    Dl[Nz-1]=-De*LL
    B[Nz] = Daq*dz*S

    # Bottom BC
    Dm[1]= 1.0
    Du[1]=-1.0
    B[1] = 0.0

    local iter
    for outer iter=1:100
        
        # Interior Points
        delta=1e-3
        Sb_p=Sb.+delta
        Sb_m=Sb.-delta
        gSb_p=g(Sb_p)
        gSb_m=g(Sb_m)
        gSb=g(Sb)
        ID=2:Nz-1
        JD=2:Nz-1
        VD=zeros(Nz-2)
        for i=2:Nz-1 
            dgds=(gSb_p[i]-gSb_m[i])/(Sb_p[i]-Sb_m[i])   
            Dm[i]=2.0 + dz^2*dgds
            B[i]=dz^2*(Sb[i]*dgds-gSb[i])
        end

        # Create tridiagonal array
        A=Tridiagonal(Dl,Dm,Du)

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
            
    return Cs,Sb,bflux,flux,iter
end