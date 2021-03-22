using SparseArrays

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
    B=zeros(Nz)

    # Create sparse matrix A 
    # Lower diagonals
    IL=2:Nz-1
    JL=1:Nz-2
    VL=-1.0*ones(Nz-2);
    # Upper diagonals
    IU=2:Nz-1
    JU=3:Nz
    VU=-1.0*ones(Nz-2);
    # Top BC
    IBt=[Nz          ,Nz    ]
    JBt=[Nz          ,Nz-1  ]
    VBt=[De*LL+Daq*dz,-De*LL]
    B[Nz]     = Daq*dz*S
    # Bottom BC
    IBb=[1  ,   1]
    JBb=[1  ,   2]
    VBb=[1.0,-1.0]
    B[1]  = 0.0

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
            VD[i-1]=2.0 + dz^2*dgds
            B[i]=dz^2*(Sb[i]*dgds-gSb[i])
        end

        # Concatenate A matrix 
        A=sparse(vcat(IL,IU,IBt,IBb,ID),vcat(JL,JU,JBt,JBb,JD),vcat(VL,VU,VBt,VBb,VD))

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