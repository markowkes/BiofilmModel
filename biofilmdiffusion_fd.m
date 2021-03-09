function [Cs,Sb,bflux,flux]=biofilmdiffusion_fd(Sbold,S,Nz,dz,param)
%% This function models the diffusion of a substrate within the biofilm
%This Function will take tank conditions (So,Xb,LL) and various growth factors (Yxs,De,Km,Daq) model the diffusion of
% substrates into the biofilm over the grid . The results of this uptake will be used to
% model the manner in which tank conditions reach equilibrium

Sb=Sbold; %preallocate array
zeroLL=1e-10; %[m] condition to consider zero thickness boundary layer
tol=1e-8; %tolerance for conversion

% Get variables out out of param
mumax=param.mumax;
Km=param.Km;
Xb=param.Xb;
Yxs=param.Yxs;
De=param.De;
LL=param.LL;
Daq=param.Daq;

%Iterations
for iter=1:100
    
    % Iteratively solves the equation in this form (and multipled by dz^2)
    %
    % Sb_{i+1} - 2 Sb_i + Sb_{i+1}     mu_max Sb_i      Xb
    % ----------------------------  = ------------ * -------
    %             dz^2                 Km + Sbold     Yxs De
    %
    A = diag((-2-dz^2*mumax./(Km+Sbold(:))*Xb/(Yxs*De)), 0) ... % Main  diagonal
        + diag( 1*ones(1,Nz-1),-1) ... % Lower diagonal
        + diag( 1*ones(1,Nz-1), 1) ;   % Upper diagonal
    B = zeros(Nz,1);
    
    %B(2:Nz+1)=mu(Sbold(2:Nz-1),param)*param.Xb/(param.Yxs*param.De);
    
    % First row - BC at bottom of biofilm
    A(1,1)=1; A(1,2)=-1; B(1,1)=0; % No flux condition
        
    % Last row - BC at top of biofilm
    if LL>zeroLL
        % Flux match condition
        A(Nz,Nz  )= De/dz+Daq/LL;
        A(Nz,Nz-1)=-De/dz;
        B(Nz,1)   =Daq/LL*S;
    else
        % Set concentration at top to S
        A(Nz,:)=0; A(Nz,Nz)=1; B(Nz,1)=S;
    end
    
    % Solve for new concentration
    Sb=A\B;
    
               
    %Non Zero Condition
    Sb(Sb < 0) = 0;
    
    % Check if converged
    if max(abs(Sb-Sbold))<tol
        break
    end

    % Transfer solution for next iteration
    Sbold=Sb;
end
Cs=Sb(end); %output Surface Concentration

%Flux Calculations
bflux=De *(Sb(end)-Sb(end-1))/dz; %"Biofilm Flux" at boundary LHS of provided flux matching equation
flux =Daq*(S      -Sb(end  ))/LL; %"Boundary Layer Flux" RHS of provided flux matching equation    
end