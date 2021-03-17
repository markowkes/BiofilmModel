using Plots
pyplot()

include("diffusion.jl")
include("cases.jl")
include("mu.jl")
param=cases(1)
#param.LL=0
Nz=50
Sbold=LinRange(0,5,Nz)
S=10
dz=1e-7
Cs,Sb,bflux,flux=diffusion(Sbold,S,Nz,dz,param)

z=LinRange(0,param.Lf,Nz)
print("Makeing Plot")
plt=plot(z,Sb)
display(plt)
