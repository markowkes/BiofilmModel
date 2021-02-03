function [S]=Diffusion(Lf,So,mewmax,Xb,Yxs,De)
%This Function will take initial tank conditions and model the diffusion of
% substrates into the biofilm. The results of this uptake will be used to
% model the manner in which tank conditions reach equilibrium

%Setup for solution
N=100; % *** Make this vary with biofilm thickness?  Or Not worth time??

%Grid of Biofilm Depth
z=linspace(0,Lf,N); %m
dz=z(2)-z(1); %m

%Substrate in Biofilm Boundary Conditions
S=zeros(1,N);
S(N)=So;

Snew=zeros(1,N);
Snew(N)=So;


%End Condition
tol=1e-6;

%Iterations

for iter=1:200
tic
        c=2:1:N-1;
        Snew(c)=(S(c+1)+S(c-1)-(mewmax*Xb*(dz^2))/(Yxs*De))/2; %Concentration of substrate at biofilm depth

        %Boundary Conditions After Iteration
        Snew(1)=S(2);
        Snew(end)=So;
        
        Snew(Snew < 0) = 0;

         S=Snew;
end
toc

% tic
% while abs(max(Snew(10)-S(10)))>tol
%         c=2:1:N-1;
%         Snew(c)=(S(c+1)+S(c-1)-(mewmax*Xb*(dz^2))/(Yxs*De))/2; %Concentration of substrate at biofilm depth
% 
%         %Boundary Conditions After Iteration
%         Snew(1)=S(2);
%         Snew(end)=So;
%         
%         Snew(Snew < 0) = 0;
%         
% end
% toc


figure(1); clf(1)
plot(z,S)
xlabel('Depth of Biofilm [m]','fontsize',20)
ylabel('Substrate Concentration within Biofilm [g/m^3]','fontsize',20)

end
