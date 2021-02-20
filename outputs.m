function outputs(t,xnew,Snew,z,bflux,Sb)
%This function organizes all the outputs into one place so that the desired
%plots can be produced in one easy step

figure(1); clf(1);
subplot(2,2,1)
plot1=plot(tnew(1:i),xnew(1:i));
xnew=xnew(1:i+1);
set(plot1,'XData',tnew,'YData',xnew)
title('Biomass Concentration For Tank')
ylabel('Concentration [g]')
xlabel('Time [s]')
hold on
subplot(2,2,2)
plot2=plot(tnew(1:i),Snew(1:i));
Snew=Snew(1:i+1);
set(plot2,'XData',tnew,'YData',Snew)
title('Substrate Concentrations For Tank')
ylabel('Concentration [g]')
xlabel('Time [s]')
% hold on
% subplot(2,2,3)
% plot3=plot(tnew(1:i),bflux(1:i));
% bflux=bflux(i+1);
% set(plot3,'XData',tnew,'YData',bflux)
% title('Flux through Boundary Layer of Biofilm')
% ylabel('Flux [g/m^2]')
% xlabel('Time [s]')
% hold on
% subplot(2,2,4)
% plot4=plot(z(1:c),Sb(1:c));
% Sb=Sb(1:c+1);
% set(plot4,'XData',z,'YData',Sb)
% title('Substrate Concentration within Biofilm')
% ylabel('Concentration [g]')
% xlabel('Depth of Biofilm [m]')

end