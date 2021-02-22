function [plots] = outputs(t,x,S,z,bflux,Sb,plots)
%This function organizes all the outputs into one place so that the desired
%plots can be produced in one easy step

% Initialize figure if start of simulation or it doesn't exist
if (all(plots==0) || ~ishandle(1)) 
    % Figure does not exist, create it
    figure(1); clf(1);
    subplot(2,2,1)
    plots(1)=plot(t,x);
    title('Biomass Concentration For Tank')
    ylabel('Concentration [g]')
    xlabel('Time [s]')
    hold on
    subplot(2,2,2)
    plots(2)=plot(t,S);
    title('Substrate Concentrations For Tank')
    ylabel('Concentration [g]')
    xlabel('Time [s]')
    drawnow
else
    % Figure exists, update data
    set(plots(1),'XData',t,'YData',x)
    set(plots(2),'XData',t,'YData',S)
    drawnow
end

% hold on
% subplot(2,2,3)
% plot3=plot(t(1:i),bflux(1:i));
% bflux=bflux(i+1);
% set(plot3,'XData',t,'YData',bflux)
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