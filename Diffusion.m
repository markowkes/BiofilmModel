function [Sb,bflux,dz]=Diffusion(Lf,LL,So,mumax,Xb,Yxs,De)
%This Function will take initial tank conditions and model the diffusion of
% substrates into the biofilm. The results of this uptake will be used to
% model the manner in which tank conditions reach equilibrium

%Setup for solution
N=100; % *** Make this vary with biofilm thickness?  Or Not worth time??

%Grid of Biofilm Depth
z=linspace(0,Lf,N); %m
dz=z(2)-z(1); %m

%Substrate in Biofilm Boundary Conditions
Sb=zeros(1,N);
Sb(end)=So; %initially assume boundary concentration = So
Snew=zeros(1,N);


%Boundary Flux Consideration
Daq=1e-4; % m^2/s (oxygen.. online)
Kl=Daq/LL; %

%Iterations
tic
for iter=1:1000

        c=2:1:N-1; %array to run concentrations through
        Snew(c)=(Sb(c+1)+Sb(c-1)-(mumax*Xb*(dz^2))/(Yxs*De))/2; %Concentration of substrate at biofilm depth
    
        %Boundary Conditions After Iteration
        Snew(1)=Snew(2);
        Snew(end)=Sb(end);
        
        %Flux Calculations
              bflux=(Snew(end)-Snew(end-1))/dz; %Biofilm Flux at boundary
              flux=(Daq*(So-Snew(end)))/(LL*De); %Boundary Layer Flux to ensure bflux=flux
              
        %Flux Matching       
        Snew(end)=So-LL*De*(bflux/Daq);
        
        %Non Zero Condition
        Snew(Snew < 0) = 0;

         Sb=Snew;
end

figure(1); clf(1)
plot(z,Sb)
xlabel('Depth of Biofilm [m]','fontsize',20)
ylabel('Substrate Concentration within Biofilm [g/m^3]','fontsize',20)

end
