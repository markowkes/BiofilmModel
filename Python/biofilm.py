#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 15:20:11 2021

@author: mowkes
"""

def tankenvironment(t,x,S,SA,Vdet,dt,Cs,Co,param):
    ## This function describes the greater tank environment and assumed that it is well mixed
    # It calls all the necessary tank geometry, flow parameters, and specific
    # parameters describing the biofilm and uses the differential equations
    # describing the substrate and biomass concentrations through time to
    # produce plots profiling these concetrations over a set period of time
    Q=param.Q
    V=param.V
    Xb=param.Xb
    Yxs=param.Yxs
    Sin=param.Sin
    Daq=param.Daq
    LL=param.LL
    
    dxdt = lambda x,t,S,Cs,Vdet: (mu(S,param)-(Q/V))*x+Vdet*SA*Xb
    dsdt = lambda x,t,S,Cs: -((mu(S,param)*x)/Yxs)+((Q*Sin)/V)-((Q*S)/V)-(SA*((Daq/LL)*(Co-Cs)))
    
    # Packing y
    y=[x, S]
    
    f = lambda t,y: np.array(
            [dxdt(y[0],t,y[1],Cs,Vdet),
             dsdt(y[0],t,y[1],Cs)]
            )
           
    tol=1e-8
    while True:
        s1 = f(t     ,y            )
        s2 = f(t+  dt/2,y+  dt/2*s1)
        s3 = f(t+3*dt/4,y+3*dt/4*s2)
        
        tnew = t + dt
        ynew = y + dt/9*(2*s1 + 3*s2 + 4*s3)
        
        s4 = f(tnew,ynew)
        
        error = dt/72*(-5*s1 + 6*s2 + 8*s3 - 9*s4)
        
        # Update timestep
        if all(abs(error) < tol/100):
            # dt is getting very small
            dt=dt*2
        elif any(abs(error) > tol):
            # dt is too big
            dt=dt/2
        else:
            # Step completed with good dt
            break

    
    # Unpacking y
    xnew=ynew[0]
    Snew=ynew[1]
        
    return tnew,xnew,Snew,dt


def biofilmDiffusion(Sbold,S,Nz,dz,param):
    
    tol=1e-12
    
    # Get variables out out of param
    Xb=param.Xb
    Yxs=param.Yxs
    De=param.De
    LL=param.LL
    Daq=param.Daq
    
    # RHS of eqn for interior points
    g = lambda S: mu(S,param)*Xb/(Yxs*De)
    
    # Initial guess
    Sb=Sbold
    
    # Preallocate A
    A=np.zeros([Nz,Nz])
    L=-np.ones(Nz-1)
    np.fill_diagonal(A[1:,:],L) # Upper diagonal
    np.fill_diagonal(A[:,1:],L) # Lower diagonal
    
    # BC at top of biofilm
    B=np.zeros([Nz,1])
    A[-1,-1]= De*LL+Daq*dz
    A[-1,-2]=-De*LL
    B[-1]     = Daq*dz*S
    
    # BC at bottom of biofilm 
    A[0,0]= 1
    A[0,1]=-1
    B[0]  = 0
    
    
    for i in range(1,100):
        
        # Interior Points
        delta=1e-3
        Sb_p=Sb[1:Nz-1]+delta
        Sb_m=Sb[1:Nz-1]-delta
        dgds=(g(Sb_p)-g(Sb_m))/(Sb_p-Sb_m)
        D=2+dz**2*dgds
        np.fill_diagonal(A[1:-1,1:-1],D)
        B[1:-1,0]=np.transpose(dz**2*(Sb[1:-1]*dgds-g(Sb[1:-1])))
     
        # Solve for new concentration
        Sb = np.linalg.solve(A,B)
        
        # Non Zero Condition
        b=np.where(Sb < 0)
        Sb[b]=0
        
        # Check if converged
        if np.max(np.abs(Sb-Sbold))<tol:
            break
    
        # Transfer solution for next iteration
        Sbold=Sb
        
    Cs=Sb[-1,0] # output Surface Concentration

    # Flux Calculations
    bflux=De *(Sb[-1,0]-Sb[-2,0])/dz # Biofilm Flux
    flux =Daq*(S       -Sb[-1,0])/LL # Boundary Layer Flux
            
    return Cs,Sb,bflux,flux


def mu(S,param):
    mumax=param.mumax
    Km=param.Km
    mu = ((mumax*S)/(Km+S))
    return mu

def lf(Sb,Lf_old,dt,dz,param):
    # This function takes the substrate concentration at a given instant in
    # time and the old biofilm thickness to computes the growth velocity as well
    # as the detachement velocity at a given instant in time. These results are
    # used to compute a new biofilm thickness Lf
        
    Kdet=param.Kdet
    
    # Growth
    Vg=0 #initial growth velocity
    for i in range(0,Sb.size-1):
        # trapezoidal integration method for new growth velocity of biofilm
        Vg=Vg+dz*0.5*(mu(Sb[i,0],param)+mu(Sb[i+1,0],param))
    
    # Detachment - New #Velocity of mass leaving biofilm into bulk liquid
    Vdet=Kdet*Lf_old**2 
    
    # Biofilm Thickness
    Lf=Lf_old+dt*(Vg-Vdet)

    return Lf,Vdet

def cases(num):
    
    # Constants
    mumax=[20, 20, 2, 20, 20, 20]
    Km   =[3, 3, 3, 3, 3, 3000]
    Yxs  =[0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
    V    =[0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
    Q    =[1, 1, 1, 50, 1, 1]
    A    =[1, 1, 1, 1, 1, 1]
    Sin  =[25, 25, 25, 25, 25, 25]
    So   =[25, 25, 25, 25, 25, 25]
    xo   =[10, 10, 10, 10, 10, 10]
    Daq  =[4.0E-5, 4.0E-5, 4.0E-5, 4.0E-5, 4.0E-5, 4.0E-5]
    De   =[1.0E-5, 1.0E-5, 1.0E-5, 1.0E-5, 1.0E-5, 1.0E-5]
    Xb   =[20000, 20000, 20000, 20000, 20000, 20000]
    Lf  =[5.0E-6, 3.0E-4, 5.0E-6, 5.0E-6, 5.0E-6, 5.0E-6]
    LL   =[1.0E-4, 1.0E-4, 1.0E-4, 1.0E-4, 1.0E-4, 1.0E-4]
    Kdet =[1900, 1900, 1900, 1900, 190000, 1900]

    # Create param
    class params:
        pass
    param = params()
    param.mumax=mumax[num]
    param.Km   =Km[num]
    param.Yxs  =Yxs[num]
    param.V    =V[num]
    param.Q    =Q[num]
    param.A    =A[num]
    param.Sin  =Sin[num]
    param.So   =So[num]
    param.xo   =xo[num]
    param.Daq  =Daq[num]
    param.De   =De[num]
    param.Xb   =Xb[num]
    param.Lf  =Lf[num]
    param.LL   =LL[num]
    param.Kdet =Kdet[num]
    return param
    
def outputs(t,x,S,z,bflux,thickness,Sb,param):
    fig, axs = plt.subplots(2,2,figsize=(10,6))
    
    ax=axs[0,0]
    ax.plot(t,x)
    ax.set_title('Tank Concentration vs. Time')
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Biomass Concentration')
    
    ax=axs[0,1]
    ax.plot(t,S)
    ax.set_title('Tank Concentration vs. Time')
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Substrate Concentration')
    
    ax=axs[1,0]
    ax.plot(z,Sb)
    ax.set_title('Substrate Concentration in Biofilm')
    ax.set_xlabel('Position within Biofilm')
    ax.set_ylabel('Substrate Concentration')    
    fig.tight_layout()
    fig.show()
    

# Main driver
import numpy as np
import matplotlib.pyplot as plt
import time
import os
os.system('cls' if os.name == 'nt' else 'clear')

start=time.time()

num=1
param=cases(num)

# Tank Parameters + Geometry
L=0.5 #[m]
W=0.5 #[m]
H=0.4 #[m]
SA=(param.V/H)+2*((param.V/L)+(param.V/W)) #tank surface area [m^2] 

Co=param.So #substrate concentration

#Create initial biofilm grid
Nz=50 #Linear GridPoints in Biofilm
z=np.linspace(0,param.Lf,Nz) #[m] Grid of Biofilm Depth
dz=z[1]-z[0] #[m]

#Initial Boundary Conditions (in Biofilm)
Sb=np.zeros(Nz)
Sb[-1]=param.So #initially assume boundary concentration = So

#Time Constraints
tFin=20 #[days]
dt=1e-2 #Interval
N=int(tFin/dt) #Number of steps
outFreq=20000 #Number of steps between plot updates.

#Preallocation
t = np.zeros(N) #Time
x = np.zeros(N) #Biomass Concentration in bulk liquid
S = np.zeros(N) #Substrate in bulk liquid
bflux=np.zeros(N) #Boundary Layer Flux of Biofilm Preallocate
flux=np.zeros(N) #Right hand side of power point equation to ensure matching flux
thickness=np.zeros(N) #Right hand side of power point equation to ensure matching flux

#Initial Conditions
t[0]=0
x[0]=param.xo
S[0]=param.So

#Initialize plots 
outIter=outFreq-1

## Time Loop
i=0
while t[i]<tFin-dt:
        
    #Update biofilm grid as biofilm grows
    z=np.linspace(0,param.Lf,Nz) #[m] Grid of Biofilm Depth
    dz=z[1]-z[0] #[m]

    #Call on "biofilmdiffusion"
    Cs,Sb,bflux[i+1],flux[i+1]=biofilmDiffusion(Sb,S[i],Nz,dz,param)
    
    #Call on "lf"
    Lf_old=param.Lf
    param.Lf,Vdet=lf(Sb,Lf_old,dt,dz,param)

    #Call on "tankenvironment"
    t[i+1],x[i+1],S[i+1],dt=tankenvironment(t[i],x[i],S[i],SA,Vdet,dt,Cs,Co,param)
    
    thickness[i+1]=param.Lf
    
    # Reallocate arrays if needed
    if t.size == i+2:
        n=t.size+1000
        t=np.resize(t,n)
        x=np.resize(x,n)
        S=np.resize(S,n)
        bflux=np.resize(bflux,n)
        flux=np.resize(flux,n)
        thickness=np.resize(thickness,n)
   
    #Call on desired plots from 'outputs'
    outIter=outIter+1
    if outIter>=outFreq:
        outputs(t[:i+1],x[:i+1],S[:i+1],z,bflux[:i+1],thickness[:i+1],Sb,param)
        outIter=0
    
    # Update iterator
    i=i+1

# Make final figures
outputs(t[:i+1],x[:i+1],S[:i+1],z,bflux[:i+1],thickness[:i+1],Sb,param)

end=time.time()
print("Elapsed time = ",end-start)

