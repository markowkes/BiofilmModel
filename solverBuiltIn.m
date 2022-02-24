
function [t,X,S,Pb,Sb,Lf]=solverBuiltIn(param)

    % Common parameters
    Nz = param.Nz;
    Nx = param.Nx;
    Ns = param.Ns;
    
    % Arrays
    X  = zeros(Nx);     % Tank particulates (concentration)
    S  = zeros(Ns);     % Tank substrates (concentration)
    Pb = zeros(Nx,Nz);  % Biofilm particulates (volume fraction)
    Sb = zeros(Nx,Nz);  % Biofilm substrates (concentration)
    
    % Initial Conditions
    X(:)    = param.xo;
    S(:)    = param.So;
    Pb(:,:) = param.phibo;
    Sb(:,:) = param.Sbo;
    Lf      = param.Lfo;

    % Package initial condition
    yo=zeros(Nx + Ns + Nx*Nz + Nx*Nz + 1, 1);
    Nvar=0;
    N=Nx;    yo(Nvar+1:Nvar+N)=X(:);    Nvar=Nvar+N; % Tank particulates
    N=Ns;    yo(Nvar+1:Nvar+N)=S(:);    Nvar=Nvar+N; % Tank substrates
    N=Nx*Nz; yo(Nvar+1:Nvar+N)=Pb(:,:); Nvar=Nvar+N; % Biofilm particulates
    N=Ns*Nz; yo(Nvar+1:Nvar+N)=Sb(:,:); Nvar=Nvar+N; % Biofilm substrates
    N=1;     yo(Nvar+1:Nvar+N)=Lf;                   % Biofilm thicknes

    % Call ODE solver
    %ops = odeset('OutputFcn',@odeprint,'AbsTol',1e-3);
    ops = odeset('OutputFcn',@myOutputFcn,'AbsTol',1e-3);
    [t,y]=ode23s(@(t,y) RHS(t,y,param),[0,param.tFin],yo,ops);

    % Extract variables
    Nvar=0;
    N=Nx;     X =y(:,Nvar+1:Nvar+N); Nvar=Nvar+N; % Tank particulates
    N=Ns;     S =y(:,Nvar+1:Nvar+N); Nvar=Nvar+N; % Tank substrates
    N=Nx*Nz;  Pb=y(:,Nvar+1:Nvar+N); Nvar=Nvar+N; % Biofilm particulates
    N=Ns*Nz;  Sb=y(:,Nvar+1:Nvar+N); Nvar=Nvar+N; % Biofilm substrates
    N=1;      Lf=y(:,Nvar+1:Nvar+N);              % Biofilm thickness

end

%% Status of ODE solver
function status=myOutputFcn(t,y,flag) %#ok<INUSL> 
    if strcmp(flag,'init') || strcmp(flag,'done')
        % do nothing
    else
        fprintf('Time = %5.5e \n',t)
    end
    status=0;
end

%% RHS of all the ODEs
function [f]=RHS(~,y,param)
    
    % Extract variables
    Nvar=0;  Nx=param.Nx; Ns=param.Ns; Nz=param.Nz;
    N=Nx;     X =y(Nvar+1:Nvar+N); Nvar=Nvar+N; % Tank particulates
    N=Ns;     S =y(Nvar+1:Nvar+N); Nvar=Nvar+N; % Tank substrates
    N=Nx*Nz;  Pb=y(Nvar+1:Nvar+N); Nvar=Nvar+N; % Biofilm particulates
    N=Ns*Nz;  Sb=y(Nvar+1:Nvar+N); Nvar=Nvar+N; % Biofilm substrates
    N=1;      Lf=y(Nvar+1:Nvar+N);              % Biofilm thickness

    % Update grid
    grid.z  = linspace(0,Lf,param.Nz);
    grid.dz = grid.z(2) - grid.z(1);
    
    % Reshape biofilm varialbes Var(Nx/Ns, Nz)
    Pb = reshape(Pb,param.Nx,param.Nz);
    Sb = reshape(Sb,param.Ns,param.Nz);

    % Compute particulate concentration from volume fractions
    Xb=zeros(param.Nx,param.Nz);
    for j=1:param.Nx
        Xb(j,:) = param.rho(j)*Pb(j,:);
    end
    
    % Compute intermediate variables
    fluxS =computeFluxS(S,Sb,param,grid);  % Flux of substrate in biofilm
    V     =computeVel  (Pb,Sb,param,grid); % Velocity of particulates
    fluxP =computeFluxP(Pb,V,param);       % Flux of particulates in biofilm
    Vdet  = param.Kdet*Lf^2;               % Detachment velocity

    
    % Compute RHS's
    f=zeros(size(y));
    Nrhs=0;
    N=Nx;    f(Nrhs+1:Nrhs+N)=dXdt (X,S,Xb,Vdet,param);      Nrhs=Nrhs+N;  % Tank particulates
    N=Ns;    f(Nrhs+1:Nrhs+N)=dSdt (X,S,param,fluxS);        Nrhs=Nrhs+N;  % Tank substrates
    N=Nx*Nz; f(Nrhs+1:Nrhs+N)=dPbdt(Pb,Sb,fluxP,param,grid); Nrhs=Nrhs+N;  % Biofilm particulates
    N=Ns*Nz; f(Nrhs+1:Nrhs+N)=dSbdt(Xb,Sb,fluxS,param,grid); Nrhs=Nrhs+N;  % Biofilm substrates
    N=1;     f(Nrhs+1:Nrhs+N)=dLfdt(V,Vdet);                               % Biofilm thickness

end

%% Fluxes of substrate due to diffusion: F=De*dSb/dz
function [fluxS]=computeFluxS(S,Sb,param,grid)
    fluxS = zeros(param.Ns,param.Nz+1); % Fluxes on faces of cells
    for i=2:param.Nz  % Interior faces
        fluxS(:,i)= param.De(:).*(Sb(:,i)-Sb(:,i-1))/grid.dz;
    end
    % Bottom boundary - no flux condition -> nothing to do
    % Top boundary - flux matching between biofilm and boundary layer 
    Sp=(param.Daq*grid.dz/2*S+param.De*param.LL*Sb(:,param.Nz)) ...
        /(param.Daq*grid.dz/2+param.De*param.LL); 
    fluxS(:,param.Nz+1) = param.Daq*(S-Sp)/param.LL;
end

%% Velocity due to growth in biofilm
function [V]=computeVel(Pb,Sb,param,grid)
    % Velocities on faces of cells
    V=zeros(1,param.Nz+1); 
    % Start with zero velocity at wall -> integrate through the biofilm
    for i=1:param.Nz
        % Start with constant velocity (no growth)
        V(i+1)=V(i);
        % Loop over particulates in this cell and add to velocity
        for j=1:param.Nx
            V(i+1)=V(i+1)+ ...
                param.mu{j}(Sb(:,i),param)*Pb(j,i)*grid.dz/param.phi_tot;
        end
    end
end

%% Fluxes of particulate due to diffusion: F=V*phi;
function [fluxP]=computeFluxP(Pb,V,param)
    % Fluxes
    fluxP = zeros(param.Nx,param.Nz+1); % Fluxes on faces of cells
    for i=2:param.Nz  % Interior faces
        fluxP(:,i)= V(i)*(Pb(:,i-1)+Pb(:,i))/2; % V*phi_face
    end
    % Bottom boundary - no flux condition -> nothing to do
    % Top boundary - use phi in top cell
    fluxP(:,param.Nz+1) = V(param.Nz+1)*(Pb(:,param.Nz));
end

%% RHS of tank particulates
function dXdt = dXdt(X,S,Xb,Vdet,param) 
    dXdt = zeros(param.Nx,1);
    for j=1:param.Nx
        dXdt(j) = param.mu{j}(S,param)*X(j) ...        % Growth
            -     param.Q*X(j)/param.V ...             % Flow out
            +     Vdet*param.A*Xb(j,end)/param.V;      % From biofilm
    end
end

%% RHS of tank substrates
function dSdt = dSdt(X,S,param,fluxS) 
    dSdt = zeros(param.Ns,1); 
    for j=1:param.Nx                                   
        dSdt(j) = param.Q.*param.Sin(j)/param.V ...       % Flow in
            -     param.Q.*      S(j)  /param.V ...       % Flow out
            -     param.A.*fluxS(j,end)/param.V ...       % Flux into biofilm
            -     param.mu{j}(S,param)*X(j)/param.Yxs(j); % Used by growth
    end
end

%% RHS of biofilm particulates 
function dPbdt = dPbdt(Pb,Sb,fluxPb,param,grid) 
    dPbdt = zeros(param.Nx,param.Nz);
    % Loop over particulates
    for j=1:param.Nx
        % Loop over cells
        for i=1:param.Nz
            dPbdt(j,i)= ...
                - (fluxPb(i+1)-fluxPb(i))/grid.dz ... % Growth velocity
                + param.mu{j}(Sb(:,i),param)*Pb(j,i); % Growth in cell
        end
    end
    % Return RHS as a column vector
    dPbdt=reshape(dPbdt,param.Nx*param.Nz,1);
end

%% RHS of biofilm substrates 
function dSbdt = dSbdt(Xb,Sb,fluxS,param,grid)
    dSbdt = zeros(param.Ns,param.Nz);
    % Loop over substrates
    for k=1:param.Ns
        % Loop over cells
        for i=1:param.Nz
            dSbdt(k,i) = (fluxS(k,i+1)-fluxS(k,i))/grid.dz; % Diffusion
            for j=1:param.Nx                                % Used by growth
                dSbdt(k,i) = dSbdt(k,i) ...
                    - param.mu{j}(Sb(:,i),param)*Xb(j,i)/param.Yxs(j);
            end
        end
    end
    % Return RHS as a column vector
    dSbdt=reshape(dSbdt,param.Ns*param.Nz,1);
end

%% RHS of biofilm thickness
function dLfdt = dLfdt(V,Vdet)
    Vfilm = V(end);    % Growth velocity at top of biofilm
    dLfdt = Vfilm - Vdet;     % Surface Velocity 
end