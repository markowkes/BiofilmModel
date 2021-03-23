using Plots, LaTeXStrings
#pyplot()
gr()
#plotlyjs()
#gr(display_type=:inline)

function outputs(t,x,S,z,bflux,Lf,Sb,param)
    p1=plot(t,x)
    xaxis!(L"\textrm{Time [days]}")
    yaxis!(L"\textrm{Biomass Con. } [g/m^3]")

    p2=plot(t,S)
    xaxis!(L"\textrm{Time [days]}")
    yaxis!(L"\textrm{Substrate Con. } [g/m^3]")

    p3=plot(t,bflux)
    xaxis!(L"\textrm{Time [days]}")
    yaxis!(L"\textrm{Flux out of Biofilm}")
    
    p4=plot(z,Sb)
    xaxis!(L"\textrm{Thickness } [\mu m]")
    yaxis!(L"\textrm{Substrate Concentration } [g/m^3]")

    p5=plot(t,Lf)
    xaxis!(L"\textrm{Time [days]}")
    yaxis!(L"\textrm{Thickness } [m]")

    N=length(t)
    dt=t[2:N]-t[1:N-1]
    #display(dt)
    p6=plot(1:(N-1),dt)
    xaxis!(L"\textrm{Iteration }")
    yaxis!(L"\textrm{Timestep } [s]")



    myplt=plot(p1,p2,p3,p4,p5,p6,layout=(2,3),
        legend=false,
        size=(800,600),
    )
    display(myplt)
    return
end