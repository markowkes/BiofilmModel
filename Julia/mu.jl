function mu(S,param)
    mumax=param.mumax
    Km=param.Km
    
    mu = ((mumax.*S)./(Km.+S))

    return mu
end