function [Sb,bflux]=biofilmdiffusion_fd(Sbold,S,Xb,dz,t,param)
%% This function models the diffusion of a substrate within the biofilm
%This Function will take tank conditions (So,Xb,LL) and various growth factors (Yxs,De,Km,Daq) model the diffusion of
% substrates into the biofilm over the grid . The results of this uptake will be used to
% model the manner in which tank conditions reach equilibrium

Ns = param.Ns;
Nz = param.Nz;
Nx = param.Nx;

Sb=Sbold;

% Get variables out of param

Yxs=param.Yxs;
De=param.De;
LL=param.LL;
Daq=param.Daq;
tol=param.dtol;

%Iterations
iter=100; %maximum iterations

L = -1;
U = -1;

% Define D = 2*kronecker + dz^2*dgds(j,i,m)
D=@(k,i,m) (2*eq(k,m)+dz^2*dgds(k,i,m,Sb,Xb,param));

% Preallocate solution array
A = zeros(Nz*Ns, Nz*Ns);
B = zeros(1,Nz*Ns);

%Iterations
for n=1:iter
    for k=1:Ns
        
        % 1 and -1 on the dia. and upper dia. at the start of each substrate (once every Nz)
        A((k-1)*Nz+1,(k-1)*Nz+1) = 1;
        A((k-1)*Nz+1,(k-1)*Nz+2) = -1;
        
        % E and F on the dia. and upper dia. at the end of each substrate (once every Nz)
        E = De(k)*LL;
        F = De(k)*LL+Daq(k)*dz;
        A(Nz*k,Nz*k) = F;
        A(Nz*k,Nz*k-1) = E;
        
        % G at the end of each substrate in B array (once every Nz)
        G = Daq(k)*dz*S(k);
        B(Nz*k) = G;

        for i=2:Nz-1
            % L and U, lower and upper dia. interior points for each
            % substrate (2:Nz-1)
            A((k-1)*Nz+i,(k-1)*Nz+i-1) = L;
            A((k-1)*Nz+i,(k-1)*Nz+i+1) = U;
            
            % R, interior points in the B matrix for each substrate
            % (2:Nz-1)
            B((k-1)*Nz+i) = R(k,i,Ns,dz,Sb,Xb,param);
            for m=1:Ns               
                
                % D, populates dia. and off dia. interior points (2:Nz-1)
                A((k-1)*Nz+i,(m-1)*Nz+i) = D(k,i,m);
                
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

% Flux = \int_0^Lf mu(S) * xB / Yxs dz = xB/Yxs * int_0^Lf mu dz
%%% Todo - Double check this %%%
bflux = zeros(Ns,1);
for k=1:Ns
    for j=1:Nx
        for i=1:length(Sb)-1
            bflux(k)=bflux(k)+dz*((Xb(j,i  )*param.mu{j}(Sb(:,i  ),param) ...
                                  +Xb(j,i+1)*param.mu{j}(Sb(:,i+1),param))/2); %trapezoidal
        end
    end
    bflux(k)=bflux(k)/Yxs(k);
end
end

function R = R(k,i,Ns,dz,Sb,Xb,param)
    R = 0;
    for m = 1:Ns
        R = R + dz^2*dgds(k,i,m,Sb,Xb,param).*Sb(m,i);
    end
    R = R - dz^2*g(k,Sb(:,i),Xb,param);
end
% Define RHS of ODE
function g = g(k,S,Xb,param)
    g = 0;
    for j = 1:param.Nx
       g = g + param.mu{j}(S,param)*Xb(j)/(param.Yxs(k)*param.De(k));
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


