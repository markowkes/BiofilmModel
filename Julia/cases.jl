
struct params
    mumax
    Km
    Yxs
    V
    Q
    A
    Sin
    So
    xo
    Daq
    De
    Xb
    Lf
    LL
    Kdet
    SA
end

function cases(n)

    # Constants
    mumax = [20 20 2 20 20 20 20];
    Km   = [3 3 3 3 3 3000 3 ];
    Yxs  = [0.5 0.5 0.5 0.5 0.5 0.5 0.5];
    V    = [0.1 0.1 0.1 0.1 0.1 0.1 0.1];
    Q    = [1 1 1 50 1 1 1];
    A    = [1 1 1 1 1 1 1];
    Sin  = [25 25 25 25 25 25 25];
    So   = [25 25 25 25 25 25 25];
    xo   = [10 10 10 10 10 10 10 ];
    Daq  = [4.0E-5 4.0E-5 4.0E-5 4.0E-5 4.0E-5 4.0E-5 4e2];
    De   = [1.0E-5 1.0E-5 1.0E-5 1.0E-5 1.0E-5 1.0E-5 1e2];
    Xb   = [20000 20000 20000 20000 20000 20000 20000];
    Lf  = [5.0E-6 3.0E-4 5.0E-6 5.0E-6 5.0E-6 5.0E-6 5e-6];
    LL   = [1.0E-4 1.0E-4 1.0E-4 1.0E-4 1.0E-4 1.0E-4 1e-4];
    Kdet = [1900 1900 1900 1900 190000 1900 1900];

    # Tank Parameters + Geometry
    L=0.5  # [m]
    W=0.5  # [m]
    H=0.4  # [m]
    SA=(V[n]/H)+2*((V[n]/L)+(V[n]/W)) # tank surface area [m^2] 


    param = params(mumax[n],Km[n],Yxs[n],V[n],Q[n],A[n],
            Sin[n],So[n],xo[n],Daq[n],De[n],Xb[n],Lf[n],LL[n],Kdet[n],SA)
    return param
end