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
%     hold on
%     subplot(2,2,3)
%     plots(3)=plot(t,bflux);
%     title('Flux through Boundary Layer of Biofilm')
%     ylabel('Flux [g/m^2]')
%     xlabel('Time [s]')
    hold on
    subplot(2,2,4)
    plots(4)=plot(z,Sb);
    title('Substrate Concentration within Biofilm')
    ylabel('Concentration [g]')
    xlabel('Depth of Biofilm [m]')
    
    drawnow
else
    % Figure exists, update data
    set(plots(1),'XData',t,'YData',x)
    set(plots(2),'XData',t,'YData',S)
%     set(plots(3),'XData',t,'YData',bflux)
    set(plots(4),'XData',z,'YData',Sb)
    drawnow
end
end