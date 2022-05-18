function [Sb,bflux]=biofilmdiffusion_fd(Sbold,S,Nz,dz,t,param)
%% This function models the diffusion of a substrate within the biofilm
%This Function will take tank conditions (So,Xb,LL) and various growth factors (Yxs,De,Km,Daq) model the diffusion of
% substrates into the biofilm over the grid . The results of this uptake will be used to
% model the manner in which tank conditions reach equilibrium

Sb=Sbold; %preallocate array
delta=1e-3;

% Get variables out of param
Xb=param.Xb;
Yxs=param.Yxs;
De=param.De;
LL=param.LL;
Daq=param.Daq;
tol=param.dtol;

%Iterations
iter=100; %maximum iterations

% Define RHS of ODE
g=@(S) mu(S,param)*Xb/(Yxs*De);

%Iterations
for i=1:iter
    
    % Interior points
    Sb_p=Sb(2:Nz-1)+delta;
    Sb_m=max(Sb(2:Nz-1)-delta,0);
    dgds=[0,(g(Sb_p)-g(Sb_m))./(Sb_p-Sb_m),0];
    A = diag(2+dz^2*dgds, 0) ... % Main  diagonal
        + diag(-1*ones(1,Nz-1),-1) ... % Lower diagonal
        + diag(-1*ones(1,Nz-1), 1) ;   % Upper diagonal
    B = (dz^2*(Sb.*dgds-g(Sb)))';
    
    % First row - No flux BC at bottom of biofilm
    A(1,1)=1; A(1,2)=-1; B(1,1)=0; 
        
    % Last row - Flux match BC at top of biofilm
    A(Nz,Nz  )= De*LL+Daq*dz;
    A(Nz,Nz-1)=-De*LL;
    B(Nz,1)   =Daq*dz*S;
    
    % Solve for new concentration
    Sb=(A\B)';
               
    % Non Zero Condition
    Sb(Sb < 0) = 0;
    
    % Check if converged
    if max(abs(Sb-Sbold))<tol
        break
    else
        if i==iter
             fprintf('Diffusion Unable to Converge at time %3.8f\n',t)
        end
    end
    
    % Transfer solution for next iteration
    Sbold=Sb;
end


% Flux = \int_0^Lf mu(S) * xB / Yxs dz = xB/Yxs * int_0^Lf mu dz
bflux = 0;
for i=1:length(Sb)-1
    bflux=bflux+dz*((mu(Sb(i),param)+mu(Sb(i+1),param))/2); %trapezoidal 
end
bflux=Xb/Yxs*bflux;

end