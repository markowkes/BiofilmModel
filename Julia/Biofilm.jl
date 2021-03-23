include("cases.jl")
include("diffusion.jl")
include("mu.jl")
include("updateLf.jl")
include("tankenvironment.jl")
include("outputs.jl")


function MainDriver()

    start = time()
    num = 7
    param = cases(num)

    # Tank Parameters + Geometry
    L = 0.5; # [m]
    W = 0.5; # [m]
    H = 0.4; # [m]
    SA = (param.V / H) + 2 * ((param.V / L) + (param.V / W)); # tank surface area [m^2] 

    Co = param.So; # substrate concentration

    # Create initial biofilm grid
    Nz = 50; # Linear GridPoints in Biofilm
    z = LinRange(0, param.Lf, Nz); # [m] Grid of Biofilm Depth

    # Initial Boundary Conditions (in Biofilm)
    Sb = zeros(Nz);
    Sb[Nz] = param.So; # initially assume boundary concentration = So

    # Time Constraints
    tFin = 30; # [days]
    dt = 1e-2; # Interval
    N = floor(Int, tFin / dt) # Number of steps
    outFreq = 2000 # Number of steps between plot updates.

    # Preallocation
    t = zeros(N); # Time
    x = zeros(N); # Biomass Concentration in bulk liquid
    S = zeros(N); # Substrate in bulk liquid
    bflux = zeros(N); # Boundary Layer Flux of Biofilm Preallocate
    flux = zeros(N); # Right hand side of power point equation to ensure matching flux
    Lf = zeros(N); # Right hand side of power point equation to ensure matching flux

    # Initial Conditions
    t[1] = 0;
    x[1] = param.xo;
    S[1] = param.So;
    Lf[1] = param.Lf

    # Initialize plots 
    outIter = outFreq - 1;
    plots = 0; titles = 0;

    ## Time Loop
    print("Starting time loop \n")
    i = 1;
    while t[i] < tFin - dt
            
        # Update biofilm grid as biofilm grows
        z = LinRange(0, Lf[i], Nz); # [m] Grid of Biofilm Depth
        dz = z[2] - z[1]; # [m]

        # Call on "biofilmdiffusion"
        Cs, Sb, bflux[i + 1], flux[i + 1],iter_diff = 
                diffusion(Sb, S[i], Nz, dz, param)
        
        # Call on "lf"
        Lf[i+1],Vdet=updateLf(Sb,Lf[i],dt,dz,param)

        # Call on "tankenvironment"
        t[i+1],x[i+1],S[i+1],dt,iter_tank=
                tankenvironment(t[i],x[i],S[i],Vdet,dt,bflux[i+1],param);
                
        # # Reallocate arrays if needed
        if length(t) == i + 1
            n = length(t) + 1000
            resize!(t, n)
            resize!(x, n)
            resize!(S, n)
            resize!(bflux, n)
            resize!(flux, n)
            resize!(Lf, n)
        end
    
        # Call on desired plots from 'outputs'
        outIter = outIter + 1;
        if outIter >= outFreq
            print("i=", i, " t=",t[i+1]," Diff Iter=",iter_diff," Tank Iter=",iter_tank,"\n")
            outputs(t[1:i+1],x[1:i+1],S[1:i+1],z,bflux[1:i+1],Lf[1:i+1],Sb,param)
            outIter = 0;
        end 

        # Update iterator
        i = i + 1

    end

    # Make final figures
    print("Making final plots \n")
    outputs(t[1:i],x[1:i],S[1:i],z,bflux[1:i],Lf[1:i],Sb,param)

    display(time()-start)
end

# Call main function
MainDriver()

