
function [t,X,S,Pb,Sb,Lf]=solverBuiltIn(param)

    % Check parameters to provide more useful error messages
    param.instantaneousDiffusion = true;
    check_param(param);

    % Common parameters
    Nz = param.Nz;
    Nx = param.Nx;
    Ns = param.Ns;
    
    % Arrays
    X  = zeros(Nx, 1);  % Tank particulates (concentration)
    S  = zeros(Ns, 1);  % Tank substrates   (concentration)
    Pb = zeros(Nx,Nz);  % Biofilm particulates (volume fraction)
    Sb = zeros(Ns,Nz);  % Biofilm substrates (concentration)
    
    % Initial Conditions
    X(:)    = param.Xo;
    S(:)    = param.So;
    for i=1:Nz
        Pb(:,i) = param.phibo;
        Sb(:,i) = param.Sbo;
    end
    Lf      = param.Lfo;

    % Package initial condition
    yo=zeros(Nx + Ns + Nx*Nz + Nx*Nz + 1, 1);
    Nvar=0;
    N=Nx;    yo(Nvar+1:Nvar+N)=X(:);    Nvar=Nvar+N; % Tank particulates
    N=Ns;    yo(Nvar+1:Nvar+N)=S(:);    Nvar=Nvar+N; % Tank substrates
    N=Nx*Nz; yo(Nvar+1:Nvar+N)=Pb(:,:); Nvar=Nvar+N; % Biofilm particulates
    if ~param.instantaneousDiffusion
        N=Ns*Nz; yo(Nvar+1:Nvar+N)=Sb(:,:); Nvar=Nvar+N; % Biofilm substrates
    end
    N=1;     yo(Nvar+1:Nvar+N)=Lf;                   % Biofilm thicknes

    % Call ODE solver
    %ops = odeset('OutputFcn',@odeprint,'AbsTol',1e-3);
    ops = odeset('OutputFcn',@myOutputFcn,'RelTol',param.tol,'AbsTol',param.tol);
    %ops = odeset('OutputFcn',@odeplot,'AbsTol',1e-4);
    %ops = odeset('OutputFcn',@odeprog,'Events',@odeabort,'AbsTol',1e-4);
    [t,y]=ode23s(@(t,y) RHS(t,y,param),[0,param.tFin],yo,ops);
    %[t,y]=ode15s(@(t,y) RHS(t,y,param),[0,param.tFin],yo,ops);
    %[t,y]=ode23t(@(t,y) RHS(t,y,param),[0,param.tFin],yo,ops);
    %[t,y]=ode23tb(@(t,y) RHS(t,y,param),[0,param.tFin],yo,ops);
    
    % Extract computed solution
    Nvar=0;
    N=Nx;     X =y(:,Nvar+1:Nvar+N); Nvar=Nvar+N; % Tank particulates
    N=Ns;     S =y(:,Nvar+1:Nvar+N); Nvar=Nvar+N; % Tank substrates
    N=Nx*Nz;  Pb=y(:,Nvar+1:Nvar+N); Nvar=Nvar+N; % Biofilm particulates
    if ~param.instantaneousDiffusion
        N=Ns*Nz;  Sb=y(:,Nvar+1:Nvar+N); Nvar=Nvar+N; % Biofilm substrates
    end
    N=1;      Lf=y(:,Nvar+1:Nvar+N);              % Biofilm thickness

    % Reshape and return last biofilm values: Var(Nx/Ns, Nz)
    Pb = reshape(Pb(end,:),Nx,Nz);

    % Substrate in biofilm
    if param.instantaneousDiffusion
        % Compute particulate concentration from volume fractions
        Xb=zeros(param.Nx,param.Nz);
        for j=1:param.Nx
            Xb(j,:) = param.rho(j)*Pb(j,:);
        end
        % Solve for final substrate concentrations in biofilm
        grid.z  = linspace(0,Lf(end),param.Nz+1);
        grid.dz = grid.z(2) - grid.z(1);
        Sb = biofilmdiffusion_fd(S(end,:),Xb,param,grid);
    else
        Sb = reshape(Sb(end,:),Ns,Nz);
    end
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
    if ~param.instantaneousDiffusion
        N=Ns*Nz;  Sb=y(Nvar+1:Nvar+N); Nvar=Nvar+N; % Biofilm substrates
    end
    N=1;      Lf=y(Nvar+1:Nvar+N);              % Biofilm thickness

    % Update grid
    grid.z  = linspace(0,Lf,param.Nz+1);
    grid.dz = grid.z(2) - grid.z(1);
    
    % Reshape biofilm varialbes Var(Nx/Ns, Nz)
    Pb = reshape(Pb,param.Nx,param.Nz);
    if ~param.instantaneousDiffusion
        Sb = reshape(Sb,param.Ns,param.Nz);
    end

    % Compute particulate concentration from volume fractions
    Xb=zeros(param.Nx,param.Nz);
    for j=1:param.Nx
        Xb(j,:) = param.rho(j)*Pb(j,:);
    end
    
    % Compute intermediate variables
    if param.instantaneousDiffusion
        Sb    = biofilmdiffusion_fd(S,Xb,param,grid); % Diffusion of substrates into biofilm
    end
    mu    = computeMu(Sb,param);        % Growthrate in biofilm
    fluxS = computeFluxS(S,Sb,param,grid);  % Flux of substrate in biofilm
    V     = computeVel  (mu,Pb,param,grid); % Velocity of particulates
    fluxP = computeFluxP(Pb,V,param);       % Flux of particulates in biofilm
    Vdet  = param.Kdet*Lf^2;                % Detachment velocity

    
    % Compute RHS's
    f=zeros(size(y));
    Nrhs=0;
    N=Nx;    f(Nrhs+1:Nrhs+N)=dXdt (X,S,Xb,Vdet,param);      Nrhs=Nrhs+N;  % Tank particulates
    N=Ns;    f(Nrhs+1:Nrhs+N)=dSdt (X,S,param,fluxS);        Nrhs=Nrhs+N;  % Tank substrates
    N=Nx*Nz; f(Nrhs+1:Nrhs+N)=dPbdt(mu,Pb,fluxP,param,grid); Nrhs=Nrhs+N;  % Biofilm particulates
    if ~param.instantaneousDiffusion
        N=Ns*Nz; f(Nrhs+1:Nrhs+N)=dSbdt(mu,Xb,fluxS,param,grid); Nrhs=Nrhs+N;  % Biofilm substrates
    end
    N=1;     f(Nrhs+1:Nrhs+N)=dLfdt(V,Vdet);                               % Biofilm thickness

end

%% Growthrate for each particulate in biofilm
function [mu]=computeMu(Sb,param)
    mu=zeros(param.Nx,param.Nz);
    % Loop over particulates
    for j=1:param.Nx
        mu(j,:)=param.mu{j}(Sb,param);
    end
end

%% Fluxes of substrate due to diffusion: F=De*dSb/dz
function [fluxS]=computeFluxS(S,Sb,param,grid)
    fluxS = zeros(param.Ns,param.Nz+1); % Fluxes on faces of cells
    for i=2:param.Nz  % Interior faces
        fluxS(:,i)= param.De(:).*(Sb(:,i)-Sb(:,i-1))/grid.dz;
    end
    % Bottom boundary - no flux condition -> nothing to do
    % Top boundary - flux matching between biofilm and boundary layer 
    S_top=(param.Daq*grid.dz/2.*S+param.De*param.LL.*Sb(:,param.Nz)) ...
        ./(param.Daq*grid.dz/2+param.De*param.LL); 
    fluxS(:,param.Nz+1) = param.De.*(S_top-Sb(:,param.Nz))/(grid.dz/2);
end

%% Velocity due to growth in biofilm
function [V]=computeVel(mu,Pb,param,grid)
    % Velocities on faces of cells
    V=zeros(1,param.Nz+1); 
    % Start with zero velocity at wall -> integrate through the biofilm
    for i=1:param.Nz
        % Add growth of particulates in this cell to velocity
        V(i+1)=V(i) + sum(mu(:,i).*Pb(:,i)*grid.dz/param.phi_tot);
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
    for k=1:param.Ns                                 
        dSdt(k) = param.Q.*param.Sin(k)/param.V ...   % Flow in
            -     param.Q.*      S(k)  /param.V ...   % Flow out
            -     param.A.*fluxS(k,end)/param.V;      % Flux into biofilm
        for j=1:param.Nx                              % Used by growth
            dSdt(k) = dSdt(k) - param.mu{j}(S,param)*X(j)/param.Yxs(j,k); 
        end
    end
end

%% RHS of biofilm particulates 
function dPbdt = dPbdt(mu,Pb,fluxPb,param,grid) 
    netFlux= (fluxPb(:,2:end)-fluxPb(:,1:end-1))/grid.dz; % Flux in/out
    growth = mu.*Pb;                                      % Growth
    dPbdt  = growth - netFlux;
    % Return RHS as a column vector
    dPbdt=reshape(dPbdt,param.Nx*param.Nz,1);
end

%% RHS of biofilm substrates 
function dSbdt = dSbdt(mu,Xb,fluxS,param,grid)
    netFlux= (fluxS(:,2:end)-fluxS(:,1:end-1))/grid.dz; % Diffusion flux
    growth = zeros(param.Ns,param.Nz);
    for k=1:param.Ns
        for j=1:param.Nx
            growth(k,:) = mu(j,:).*Xb(j,:)./param.Yxs(j,k); % Used by growth
        end
    end
    dSbdt = netFlux - growth; 
    % Return RHS as a column vector
    dSbdt=reshape(dSbdt,param.Ns*param.Nz,1);
end

%% RHS of biofilm thickness
function dLfdt = dLfdt(V,Vdet)
    Vfilm = V(end);    % Growth velocity at top of biofilm
    dLfdt = Vfilm - Vdet;     % Surface Velocity 
end

%% Check parameters
function check_param(param)
    Nx=param.Nx;
    Ns=param.Ns;
    Nz=param.Nz;

    if param.tFin<=0
        error('param.tFin should be a positive amount of time')
    end

    if param.V<=0
        error('param.V should be a positive volume')
    end

    if param.A<=0
        error('param.A should be a positive area')
    end

    if param.Q<0
        error('param.Q should be a non-negative flow rate')
    end

    if length(param.Xo)~=Nx
        error(['param.Xo should have Nx=',num2str(Nx),' particulate tank ICs'])
    end

    if length(param.So)~=Ns
        error(['param.So should have Ns=',num2str(Ns),' substrate tank ICs'])
    end

    if length(param.Sin)~=Ns
        error(['param.Sin should have Ns=',num2str(Ns),' substrate concentrations'])
    end

    if length(param.LL)<0
        error('param.LL should be a non-negative boundary layer thickness')
    end

    if param.Nz<2 || param.Nz ~= floor(param.Nz)
        error('param.Nz should be an integer >= 2')
    end

    if length(param.phibo)~=Nx
        error(['param.phibo should have Nx=',num2str(Nx),' particulate ICs for biofilm'])
    end

    if length(param.Sbo)~=Ns
        error(['param.Sbo should have Ns=',num2str(Ns),' substrate ICs for biofilm'])
    end

    if length(param.Lfo)~=1 || param.Lfo<1e-12
        error('param.Lfo should have one number that is larger than 1e-12 (an already very small number)')
    end

    if size(param.Yxs,1)~=Nx || size(param.Yxs,2)~=Ns
        error(['param.Yxs should be of size', Nx,' x ',Ns])
    end
    if any(param.Yxs<eps) 
        error('param.Yxs should have all non-zero entries')
    end

    if length(param.Daq) ~= Ns
        error(['param.Daq should have Ns=',num2str(Ns),' diffusion coefficients'])
    end

    if length(param.De) ~= Ns
        error(['param.De should have Ns=',num2str(Ns),' diffusion coefficients'])
    end

    if length(param.rho) ~= Nx
        error(['param.rho should have Nx=',num2str(Nx),' densities'])
    end

    if length(param.Kdet) ~= 1
        error('param.Kdet should have 1 detachment coefficient')
    end

    if length(param.mu) ~= Nx
        error(['param.mu should have Nx=',num2str(Ns),' growthrate equations'])
    end

    Stest=rand(Ns,Nz);
    for j=1:Nx
        if ~isequal(size(param.mu{j}(Stest,param)),[1 Nz])
            error(['mu{',num2str(j),'}(Sb,param) returns a matrix of size '...
                ,num2str(size(param.mu{j}(Stest,param))),[', ' ...
                'it should return a matrix of size 1 x '],num2str(Nz)])
        end
    end


end