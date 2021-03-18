using Plots
#pyplot()
#gr()
plotlyjs()

function outputs(t,x,S,z,bflux,Lf,Sb,param)
    p1=plot(t,x)
    xaxis!("Time")
    yaxis!("Biomass")

    p2=plot(t,S)
    xaxis!("Time")
    yaxis!("Substrate Concentration")

    p3=plot(z,Sb)
    xaxis!("Thickness")
    yaxis!("Substrate Concentration")

    p4=plot(t,Lf)
    xaxis!("Time")
    yaxis!("Thickness")

    plt=plot(p1,p2,p3,p4,layout=(2,2))
    display(plt)
    return
end