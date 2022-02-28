function plotSolution(t,X,S,Pb,Sb,Lf,param)

%% Time dependent outputs 

% Tank particulates
figure(1); clf(1)
subplot(2,3,1)
plot(t,X)
xlabel('Time')
ylabel('Tank Particulate Concentrations')
set(gca,'Fontsize',16)

% Tank substrates
subplot(2,3,2)
plot(t,S)
xlabel('Time')
ylabel('Tank Substrate Concentrations')
set(gca,'Fontsize',16)

% Biofilm thickness
subplot(2,3,3)
plot(t,Lf*1e6)
xlabel('Time')
ylabel('Biofilm Thickness (um)')
set(gca,'Fontsize',16)

%% Variables in biofilm (vary spatially at last time)
Nz=param.Nz;
z  = linspace(0,Lf(end),Nz+1); % Grid points
zm = 0.5*(z(1:Nz)+z(2:Nz+1));  % Grid cells
dz = z(2)-z(1);                % Grid spacing

% Biofilm particulates
subplot(2,3,4)
hold on
for j=1:param.Nx
    plot(zm,Pb(j,:))
end
xlabel('Location in Biofilm')
ylabel('Biofilm Particulate Vol. Fraction')
ylim([0,1])
set(gca,'Fontsize',16)

% Biofim substrates
subplot(2,3,5)
hold on
% Add point at the top of biofilm (flux matching condition)
S_top=(param.Daq*dz/2.*S(end,:)'+param.De*param.LL.*Sb(:,param.Nz)) ...
        ./(param.Daq*dz/2+param.De*param.LL); 
for k=1:param.Ns
    plot([zm,z(end)],[Sb(k,:),S_top(k)])
end
xlabel('Location in Biofilm')
ylabel('Biofilm Substrate Concentrations')
set(gca,'Fontsize',16)
