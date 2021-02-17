clear;
clc;

%Steady State Test Conditions
mumax=20; %d-1	maximum specific growth rate
Km=3; %g m-3	Monod half saturation coefficient
Yxs=0.5; %gx gs-1	biomass yield coefficient on substrate
De=5.00E-05; %m2 d-1	effective diffusion coefficient of substrate in biofilm
Xb=20000; %g m-3	biomass density in biofilm
Lf=4.00E-4; %m	biofilm thickness
So=25; %g m-3	bulk fluid substrate concentration

% So=input('initial concentration\n');
% Lf=input('inital thickness\n');
% mewmax=input('specific growth rate\n');
% Yxs=input('biomass yield coeff on substrate\n');
% De=input('effective diffusion\n');
% Xb=input('biomass density\n');

[Sb,bflux]=Diffusion(Lf,So,mumax,Xb,Yxs,De)

[mu] = mu_function(mumax,Km,Sb)
