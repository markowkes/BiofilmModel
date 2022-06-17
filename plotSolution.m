function plotSolution(t,X,S,Pb,Sb,Lf,param,type)

if ~exist('type','var')
    type='final';
end

% Colors
switch type
    case 'update'
        Cline = {'k.','b.','r.','g.','y.'};
    case 'final'
        Cline = {'k-','b-','r-','g-','y-'};
    otherwise 
        error('Unknown type in plotSolution')
end


C = {'k','b','r','g','y'};
    
% Tank particulates
figure(1)
subplot(2,3,1)
hold on
for i=1:param.Nx
    plot(t,X(i,:),Cline{i})
end
xlabel('Time')
ylabel('Tank Particulate Concentrations')
legend(param.XNames)
%ylim([min(X,[],'all')-0.3*max(0.1,(max(X,[],'all')-min(X,[],'all'))) ...
%      max(X,[],'all')+0.3*max(0.1,(max(X,[],'all')-min(X,[],'all')))])
set(gca,'Fontsize',16)

% Tank substrates
subplot(2,3,2)
hold on
for j=1:param.Ns
    plot(t,S(j,:),Cline{j})
end
xlabel('Time')
ylabel('Tank Substrate Concentrations')
legend(param.SNames)
ylim([min(S,[],'all')-0.3*max(0.1,(max(S,[],'all')-min(S,[],'all'))) ...
    max(S,[],'all')+0.3*max(0.1,(max(S,[],'all')-min(S,[],'all')))])
set(gca,'Fontsize',16)

% Biofilm thickness
subplot(2,3,3)
hold on
plot(t,Lf*1e6,Cline{1})
xlabel('Time')
ylabel('Biofilm Thickness (um)')
%ylim(1e6*[min(Lf,[],'all')-0.3*max(0.1,(max(Lf,[],'all')-min(Lf,[],'all'))) ...
%          max(Lf,[],'all')+0.3*max(0.1,(max(Lf,[],'all')-min(Lf,[],'all')))])
set(gca,'Fontsize',16)

%% Variables in biofilm (vary spatially at last time)
Nz=param.Nz;
z  = linspace(0,Lf(end),Nz+1); % Grid points
zm = 0.5*(z(1:Nz)+z(2:Nz+1));  % Grid cells
dz = z(2)-z(1);                % Grid spacing

% Biofilm particulates
subplot(2,3,4); cla()
hold on
for j=1:param.Nx
    plot(zm,Pb(j,:),C{j})
end
xlabel('Location in Biofilm')
ylabel('Biofilm Particulate Vol. Fraction')
legend(param.XNames)
%ylim([min(min(Pb))-0.3*(max((max(Pb))-min(min(Pb)))) ...
%    max(max(Pb))+0.3*(max(max(Pb))-min(min(Pb)))])
set(gca,'Fontsize',16)

% Biofim substrates
subplot(2,3,5); cla()
hold on
% Add point at the top of biofilm (flux matching condition)
S_top=(param.Daq*dz/2.*S(:,end)'+param.De*param.LL.*Sb(:,param.Nz)) ...
    ./(param.Daq*dz/2+param.De*param.LL);
for k=1:param.Ns
    plot([zm,z(end)],[Sb(k,:),S_top(k)],C{k})
end
xlabel('Location in Biofilm')
ylabel('Biofilm Substrate Concentrations')
legend(param.SNames)
ylim([min(min(Sb))-0.3*(max((max(Sb))-min(min(Sb)))) ...
    max(max(Sb))+0.3*(max(max(Sb))-min(min(Sb)))])
set(gca,'Fontsize',16)


