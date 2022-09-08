function [Sb,Sflux]=biofilmdiffusion_fd(t,S,Xb,Lf,param,grid)
%% This function models the diffusion of a substrate within the biofilm
%This Function will take tank conditions (So,Xb,LL) and various growth factors (Yxs,De,Km,Daq) model the diffusion of
% substrates into the biofilm over the grid . The results of this uptake will be used to
% model the manner in which tank conditions reach equilibrium

Ns = param.Ns;
Nz = param.Nz;

Sbold=zeros(Ns,Nz);
Sb =S.*ones(Ns,Nz);

% Get variables out of param
De=param.De;
LL=param.LL;
Daq=param.Daq;
dz = grid.dz;
tol=param.tol;

%Iterations
iter=100; %maximum iterations

L = -1;
U = -1;

%Iterations
for n=1:iter

    % Preallocate solution array
    A = zeros(Nz*Ns, Nz*Ns);
    B = zeros(Nz*Ns,1);

    % Precompute g(k) at z(i)  %%%% Need to vectoriz %%%%
    g=zeros(Nz,Ns);
    for k=1:Ns
        g(:,k)=compute_g(t,k,Sb,Xb,Lf,param,grid);
    end

    % Precompute dg(k)/ds(m) at z(i)
    dgds=zeros(Ns,Nz,Ns);
    for k=1:Ns
        for m=1:Ns
            dgds(m,:,k)=compute_dgds(t,k,m,Sb,Xb,Lf,g,param,grid);
        end
    end
    

    % Loop over substrates
    for k=1:Ns
        
        % Bottom of biofilm
        i=1; % Index
        d=(k-1)*Nz+i; % Row for this substrate 
        A(d,d  ) = L + 2;  % Add lower diagonal to main diagonal (Neuman)
        A(d,d+1) = U; 
        B(d,1  ) = R(k,i,Ns,dz,Sb,g,dgds); % RHS
        
        % Top of biofilm
        i=Nz; % Index
        d=(k-1)*Nz+i; % Row for this substrate 
        % Top of biofilm - flux matching condition between biofilm and tank
        A(d,d  ) = ((3*Daq(k)*dz + 2*De(k)*LL))/(Daq(k)*dz + 2*De(k)*LL);
        A(d,d-1) = L; 
        B(Nz*k) = R(k,i,Ns,dz,Sb,g,dgds) + (2*Daq(k)*S(k)*dz)/(Daq(k)*dz + 2*De(k)*LL);

        % Interior points
        for i=2:Nz-1
            d=(k-1)*Nz+i; % Row for this substrate 
            A(d,d-1) = L;
            A(d,d+1) = U;
            A(d,d  ) = 2;
            B(d,1  ) = R(k,i,Ns,dz,Sb,g,dgds);
        end

        % All points
        for i=1:Nz
            d=(k-1)*Nz+i; % Row for this substrate 
            for m=1:Ns               
                % D, populates dia. and off dia. interior points
                A(d,(m-1)*Nz+i) = A(d,(m-1)*Nz+i) ...
                    + dz^2*dgds(m,i,k);
            end
        end
    end

    % Solve for new concentration
    Sb=(A\B)';

    % Non Zero Condition
    Sb(Sb < 0) = 0;

    % Reshape from vector to a more readable matrix
    Sb = reshape(Sb,[Nz,Ns]);
    Sb = Sb';
    
    % Check if converged
    if max(abs(Sb-Sbold))<tol
        break
    else
        if n==iter
             fprintf('Diffusion Unable to Converge at time %3.8f\n',t)
        end
    end
    % Save current iteration Sb value
    Sbold = Sb;
end

% Compute Sflux at top of biofilm
S_top = (Daq*dz/2.*S + De*LL.*Sb(:,end))./(Daq*dz/2 + De*LL);
Sflux = De.*(S_top - Sb(:,end))/(dz/2);
end

% k -> substrate
% i -> location in biofilm
function R = R(k,i,Ns,dz,Sb,g,dgds)
    R = 0;
    for m = 1:Ns
        R = R + dz^2*dgds(m,i,k).*Sb(m,i);
    end
    R = R - dz^2*g(i,k);
end
% Define RHS of ODE
function g = compute_g(t,k,Sb,Xb,Lf,param,grid)
    theavi = mod(t, 1);
    g = sum(param.mu(Sb,Xb,Lf,theavi,grid.z(1:end-1),param).*Xb./(param.Yxs(:,k)*param.De(k)),1);
end
% Define dgds = (g(Sb+)-g(Sb-))/dS
function dgds = compute_dgds(t,k,m,Sb,Xb,Lf,g,param,grid) 
    % Define Sb plus and Sb minus, delta is added/subtracted to Sb(i,m)
    delta=1e-3; 
    Sb_p = Sb; Sb_p(m,:)=Sb_p(m,:)+delta;
    % Compute g at plus and minus points
    gp=compute_g(t,k,Sb_p,Xb,Lf,param,grid);
    gm=g(:,k)';
    % Compute derivative
    dgds=(gp-gm)/delta;
end


