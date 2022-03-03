function [Sb,Sflux]=biofilmdiffusion_fd(S,Xb,param,grid)
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
D=@(k,i,m,Sb,Xb) (dz^2*dgds(k,i,m,Sb,Xb,param));

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
        B(d,1  ) = R(k,i,Ns,dz,Sb,Xb,param); % RHS
        
        % Top of biofilm
        i=Nz; % Index
        d=(k-1)*Nz+i; % Row for this substrate 
        % Top of biofilm - flux matching condition between biofilm and tank
        A(d,d  ) = ((3*Daq(k)*dz + 2*De(k)*LL))/(Daq(k)*dz + 2*De(k)*LL);
        A(d,d-1) = L; 
        B(Nz*k) = R(k,i,Ns,dz,Sb,Xb,param) + (2*Daq(k)*S(k)*dz)/(Daq(k)*dz + 2*De(k)*LL);

        % Interior points
        for i=2:Nz-1
            d=(k-1)*Nz+i; % Row for this substrate 
            A(d,d-1) = L;
            A(d,d+1) = U;
            A(d,d  ) = 2;
            B(d,1  ) = R(k,i,Ns,dz,Sb,Xb,param);
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
function R = R(k,i,Ns,dz,Sb,Xb,param)
    R = 0;
    for m = 1:Ns
        R = R + dz^2*dgds(k,i,m,Sb,Xb,param).*Sb(m,i);
    end
    R = R - dz^2*g(k,Sb(:,i),Xb(:,i),param);
end
% Define RHS of ODE
function g = g(k,S,Xb,param)
    g = 0;
    for j = 1:param.Nx
       g = g + param.mu{j}(S,param)*Xb(j)/(param.Yxs(j,k)*param.De(k));
    end
end
% Define dgds = (g(Sb+)-g(Sb-))/dS
function dgds = dgds(k,i,m,Sb,Xb,param) 
    % Define Sb plus and Sb minus, delta is added/subtracted to Sb(i,m)
    delta=1e-3; 
    Sb_p =        Sb(:,i) + delta.*transpose(eq(1:param.Ns,m)) ;
    Sb_m = max(0, Sb(:,i) - delta.*transpose(eq(1:param.Ns,m)));
    % Compute g at plus and minus points
    gp=g(k,Sb_p,Xb(:,i),param);
    gm=g(k,Sb_m,Xb(:,i),param);
    % Compute derivative
    dgds=(gp-gm)/max(Sb_p-Sb_m);
end


