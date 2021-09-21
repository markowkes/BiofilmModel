function [Sb,bflux]=biofilmdiffusion_fd(Sbold,S,Nz,dz,t,param)
%% This function models the diffusion of a substrate within the biofilm
%This Function will take tank conditions (So,Xb,LL) and various growth factors (Yxs,De,Km,Daq) model the diffusion of
% substrates into the biofilm over the grid . The results of this uptake will be used to
% model the manner in which tank conditions reach equilibrium

Ns = size(S,1);

Sb=Sbold;
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

L = -1;
U = -1;
j=1;
% Define RHS of ODE
g=@(k,S) mu(j,S,param)*Xb/(Yxs*De(k));

% Define Sb plus and Sb minus, delta is added/subtracted to Sb(i,m)
Sb_p=@(m,i) Sbold(m,:) + delta.*transpose(eq(1:Ns,m));
Sb_m=@(m,i) max(0, Sbold(m,:) - delta.*transpose(eq(1:Ns,m)));

% Define dgds = (g(Sb+)-g(Sb-))/dS
dgds=@(k,i,m) (g(k,Sb_p(m,i))-g(k,Sb_m(m,i)))/max(max(Sb_p(m,i) - Sb_m(m,i)));

% Define D = 2*kronecker + dz^2*dgds(j,i,m)
D=@(k,i,m) (2*eq(k,m)+dz^2*dgds(k,i,m));

% Preallocate solution array
A = zeros(Nz*Ns, Nz*Ns);

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
            B((k-1)*Nz+i) = R(k,i,Ns,dz,Sbold,dgds,g,Sb_p,Sb_m,delta);
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
    Sb = reshape(Sb,[Ns,Nz]);
    
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
    bflux=bflux+dz*((mu(1,Sb(:,i),param)+mu(1,Sb(:,i+1),param))/2); %trapezoidal 
end
bflux=Xb/Yxs.*bflux;




function r = R(k,i,Ns,dz,Sbold,dgds,g,Sb_p,Sb_m,delta)
    r = 0;
    for m = 1:Ns
        r = r + dz^2*dgds(k,i,m).*Sbold(m,i);
    end
    r = r - dz^2*g(k,Sbold(:,i));
