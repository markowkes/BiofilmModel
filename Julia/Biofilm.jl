include("cases.jl")
include("diffusion.jl")
include("mu.jl")
include("Lf.jl")
include("tankenvironment.jl")
include("outputs.jl")


function MainDriver()

    num = 1
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
    tFin = 20; # [days]
    dt = 1e-2; # Interval
    N = floor(Int, tFin / dt) # Number of steps
    outFreq = 20000 # Number of steps between plot updates.

    # Preallocation
    t = zeros(N); # Time
    x = zeros(N); # Biomass Concentration in bulk liquid
    S = zeros(N); # Substrate in bulk liquid
    bflux = zeros(N); # Boundary Layer Flux of Biofilm Preallocate
    flux = zeros(N); # Right hand side of power point equation to ensure matching flux
    thickness = zeros(N); # Right hand side of power point equation to ensure matching flux

    # Initial Conditions
    t[1] = 0;
    x[1] = param.xo;
    S[1] = param.So;

    # Initialize plots 
    outIter = outFreq - 1;
    plots = 0; titles = 0;

    ## Time Loop
    print("Starting time loop \n")
    i = 1;
    while t[i] < tFin - dt
            
        # Update biofilm grid as biofilm grows
        z = LinRange(0, param.Lf, Nz); # [m] Grid of Biofilm Depth
        dz = z[2] - z[1]; # [m]

        # Call on "biofilmdiffusion"
        Cs, Sb, bflux[i + 1], flux[i + 1] = diffusion(Sb, S[i], Nz, dz, param)
        
        # Call on "lf"
        Lf_old = param.Lf;
        param.Lf,Vdet=Lf(Sb,Lf_old,dt,dz,param)

        # Call on "tankenvironment"
        t[i+1],x[i+1],S[i+1],dt=tankenvironment(t[i],x[i],S[i],SA,Vdet,dt,Cs,Co,param);
        
        thickness[i + 1] = param.Lf;
        
        # # Reallocate arrays if needed
        if length(t) == i + 1
            n = length(t) + 1000
            resize!(t, n)
            resize!(x, n)
            resize!(S, n)
            resize!(bflux, n)
            resize!(flux, n)
            resize!(thickness, n)
        end
    
        # Call on desired plots from 'outputs'
        outIter = outIter + 1;
        if outIter >= outFreq
            print("i=", i, "\n")
            outputs(t[1:i+1],x[1:i+1],S[1:i+1],z,bflux[1:i+1],thickness[1:i+1],Sb,param)
            outIter = 0;
        end 

        # Update iterator
        i = i + 1

    end

    # Make final figures
    print("Making final plots")
    outputs(t[1:i+1],x[1:i+1],S[1:i+1],z,bflux[1:i+1],thickness[1:i+1],Sb,param)

end

# Call main function
MainDriver()

