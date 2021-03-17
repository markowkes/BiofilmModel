using Plots
pyplot()

function outputs(t,x,S,z,bflux,thickness,Sb,param)
    plt=plot(x,t)
    display(plt)
    return
end