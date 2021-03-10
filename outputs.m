function [plots,titles] = outputs(t,x,S,z,bflux,thickness,Sb,param,plots,titles)
%This function organizes all the outputs into one place so that the desired
%plots can be produced in one easy step

% Initialize figure if start of simulation or it doesn't exist
if (all(plots(:)==0) || ~ishandle(1)) 
    % Figure does not exist, create it
    figure(1); clf(1);
    subplot(2,3,1)
    plots(1,1)=plot(t,x);
    title('Biomass Concentration For Tank')
    ylabel('Concentration [g/m^3]')
    xlabel('Time [days]')
    %xlim([0,2])
    %ylim([0,25])
    hold on
    
    subplot(2,3,2)
    plots(2,1)=plot(t,S);
    title('Substrate Concentrations For Tank')
    ylabel('Concentration [g/m^3]')
    xlabel('Time [days]')
    %xlim([0,2])
    %ylim([0,25])
    hold on
    
    subplot(2,3,3)
    plots(3,1)=plot(t,bflux);
    title('Flux through Boundary Layer of Biofilm')
    ylabel('Flux [g/m^2]')
    xlabel('Time [days]')
    %xlim([0,2])
    hold on
    
    subplot(2,3,4)
    plots(4,1)=plot(z,Sb);
    hold on
    plots(4,2)=plot([z(end),z(end)+param.LL],[Sb(end),S(end)]);
    title('Substrate Concentration within Biofilm')
    ylabel('Concentration [g/m^3]')
    xlabel('Depth of Biofilm [m]')
    %xlim([0,Lf])
    
    subplot(2,3,5)
    plots(5,1)=plot(t,bflux);
    title('Thickness of biofilm vs Time')
    ylabel('Thickness [m]')
    xlabel('Time [days]')
    %xlim([0,2])
    hold on
    
    % Title
    str=sprintf('Time = %5.2f days',t(end));
    titles=sgtitle(str);
    drawnow
else
    % Figure exists, update data
    set(plots(1,1),'XData',t,'YData',x)
    set(plots(2,1),'XData',t,'YData',S)
    set(plots(3,1),'XData',t,'YData',bflux)
    set(plots(4,1),'XData',z,'YData',Sb)
    set(plots(4,2),'XData',[z(end),z(end)+param.LL],'YData',[Sb(end),S(end)]);
    set(plots(5,1),'XData',t,'YData',thickness)
    str=sprintf('Time = %5.2f days',t(end));
    set(titles,'String',str)
    drawnow
end
end