function [Sb,Sflux]=biofilmdiffusion_fd(t,S,Xb,param,grid)
%% This function models the diffusion of a substrate within the biofilm
%This Function will take tank conditions (So,Xb,LL) and various growth factors (Yxs,De,Km,Daq) model the diffusion of
% substrates into the biofilm over the grid . The results of this uptake will be used to
% model the manner in which tank conditions reach equilibrium

Ns = param.Ns;
Nz = param.Nz;

Sbold=zeros(Ns,Nz);
Sb   =zeros(Ns,Nz);

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

% Define D = 2*kronecker + dz^2*dgds(j,i,m)
D=@(k,i,m,Sb,Xb) (dz^2*dgds(t,k,i,m,Sb,Xb,param,grid));

%Iterations
for n=1:iter
    % Preallocate solution array
    A = zeros(Nz*Ns, Nz*Ns);
    B = zeros(Nz*Ns,1);

    for k=1:Ns
        
        % Bottom of biofilm
        i=1; % Index
        d=(k-1)*Nz+i; % Row for this substrate 
        A(d,d  ) = L + 2;  % Add lower diagonal to main diagonal (Neuman)
        A(d,d+1) = U; 
        B(d,1  ) = R(t,k,i,Ns,dz,Sb,Xb,param,grid); % RHS
        
        % Top of biofilm
        i=Nz; % Index
        d=(k-1)*Nz+i; % Row for this substrate 
        % Top of biofilm - flux matching condition between biofilm and tank
        A(d,d  ) = ((3*Daq(k)*dz + 2*De(k)*LL))/(Daq(k)*dz + 2*De(k)*LL);
        A(d,d-1) = L; 
        B(Nz*k) = R(t,k,i,Ns,dz,Sb,Xb,param,grid) + (2*Daq(k)*S(k)*dz)/(Daq(k)*dz + 2*De(k)*LL);

        % Interior points
        for i=2:Nz-1
            d=(k-1)*Nz+i; % Row for this substrate 
            A(d,d-1) = L;
            A(d,d+1) = U;
            A(d,d  ) = 2;
            B(d,1  ) = R(t,k,i,Ns,dz,Sb,Xb,param,grid);
        end

        % All points
        for i=1:Nz
            d=(k-1)*Nz+i; % Row for this substrate 
            for m=1:Ns               
                % D, populates dia. and off dia. interior points (2:Nz-1)
                A(d,(m-1)*Nz+i) = A(d,(m-1)*Nz+i) + D(k,i,m,Sb,Xb);
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
        if i==iter
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
function R = R(t,k,i,Ns,dz,Sb,Xb,param,grid)
    R = 0;
    for m = 1:Ns
        R = R + dz^2*dgds(t,k,i,m,Sb,Xb,param,grid).*Sb(m,i);
    end
    R = R - dz^2*g(t,k,i,Sb(:,i),Xb(:,i),param,grid);
end
% Define RHS of ODE
function g = g(t,k,i,Sb,Xb,param,grid)
    g = 0;
    theavi = mod(t, 1);
    for j = 1:param.Nx
       g = g + param.mu{j}(Sb,Xb,theavi,grid.z(i),param)*Xb(j)/(param.Yxs(j,k)*param.De(k));
    end
end
% Define dgds = (g(Sb+)-g(Sb-))/dS
function dgds = dgds(t,k,i,m,Sb,Xb,param,grid) 
    % Define Sb plus and Sb minus, delta is added/subtracted to Sb(i,m)
    delta=1e-3; 
    Sb_p =        Sb(:,i) + delta.*transpose(eq(1:param.Ns,m)) ;
    Sb_m = max(0, Sb(:,i) - delta.*transpose(eq(1:param.Ns,m)));
    % Compute g at plus and minus points
    gp=g(t,k,i,Sb_p,Xb(:,i),param,grid);
    gm=g(t,k,i,Sb_m,Xb(:,i),param,grid);
    % Compute derivative
    dgds=(gp-gm)/max(Sb_p-Sb_m);
end


