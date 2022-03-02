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
D=@(k,i,m,Sb,Xb) (2*eq(k,m)+dz^2*dgds(k,i,m,Sb,Xb,param));

% Preallocate solution array
A = zeros(Nz*Ns, Nz*Ns);
B = zeros(1,Nz*Ns);

%Iterations
for n=1:iter
    for k=1:Ns
        
        % E and F on the dia. and upper dia. at the end of each substrate (once every Nz)
%         E = -De(k)*LL;
%         F =  De(k)*LL+Daq(k)*dz;
%         A(Nz*k,Nz*k) = F;
%         A(Nz*k,Nz*k-1) = E;
        A(Nz*k,Nz*k) = ((3*Daq(k)*dz + 2*De(k)*LL))/(Daq(k)*dz + 2*De(k)*LL); % Si term
        A(Nz*k,Nz*k-1) = -1;  % Si-1 term
        
        % G at the end of each substrate in B array (once every Nz)
        %G = Daq(k)*dz*S(k);
        %B(Nz*k) = G;
        B(Nz*k) = R(k,Nz,Ns,dz,Sb,Xb,param) + (2*Daq(k)*S(k)*dz)/(Daq(k)*dz + 2*De(k)*LL);


        for i=1:Nz-1
            % L and U, lower and upper dia. interior points for each
            % substrate (2:Nz-1)
            if i~=1
                A((k-1)*Nz+i,(k-1)*Nz+i-1) = L;
            else
                A((k-1)*Nz+i,(k-1)*Nz+i  ) = L;  % add to main diagonal (BC)
            end
            A((k-1)*Nz+i,(k-1)*Nz+i+1) = U;
            
            % R, interior points in the B matrix for each substrate
            % (2:Nz-1)
            B((k-1)*Nz+i) = R(k,i,Ns,dz,Sb,Xb,param);
            for m=1:Ns               
                
                % D, populates dia. and off dia. interior points (2:Nz-1)
                A((k-1)*Nz+i,(m-1)*Nz+i) = A((k-1)*Nz+i,(m-1)*Nz+i) + D(k,i,m,Sb,Xb);
                
%                 A((j-1)*Nz+i,(m-1)*Nz+i+1) = U;
%                 A((j-1)*Nz+i,(m-1)*Nz+i-1) = L;     
            end

        end
    end

    % Solve for new concentration
    Sb=(A\B')';

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
Stop = (Daq*dz/2.*S + De*LL.*Sb(:,end))./(Daq*dz/2 + De*LL);
Sflux = (Stop - Sb(:,end))/(dz/2);
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


