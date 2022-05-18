function [plots,titles] = outputs(t,x,S,z,bflux,thickness,Sb,param,plots,titles)
%This function organizes all the outputs into one place so that the desired
%plots can be produced in one easy step

% Initialize figure if start of simulation or it doesn't exist
if (all(plots(:)==0) || ~ishandle(1)) 
    % Figure does not exist, create it
    figure(1); clf(1);
    subplot(2,3,1)
    for i=1:param.Nx
        plots(1,i)=plot(t,x(i,:));
        hold on
    end
    title('Biomass Concentration For Tank')
    ylabel('Concentration [g/m^3]')
    xlabel('Time [days]')
    hold on
    
    subplot(2,3,2)
    for i=1:param.Ns
        plots(2,i)=plot(t,S(1,:));
        hold on
    end
    title('Substrate Concentrations For Tank')
    ylabel('Concentration [g/m^3]')
    xlabel('Time [days]')
    hold on
    
    subplot(2,3,3)
    for i=1:param.Ns
        plots(3,i)=plot(t,bflux(i,:));
        hold on
    end
    title('Flux through Boundary Layer of Biofilm')
    ylabel('Flux [g/m^2]')
    xlabel('Time [days]')
    hold on
    
    subplot(2,3,4)
    for i=1:param.Ns
        plots(4,i)=plot(z,Sb(i,:));
        hold on
        %plots(4,2)=plot([z(end),z(end)+param.LL],[Sb(end),S(end)]);
    end
    title('Substrate Concentration within Biofilm')
    ylabel('Concentration [g/m^3]')
    xlabel('Depth of Biofilm [m]')
   
    
    subplot(2,3,5)
    plots(5,1)=plot(t,1E6*thickness);
    title('Thickness of biofilm vs Time')
    ylabel('Thickness [micrometers]')
    xlabel('Time [days]')
   
    hold on
    
    % Title
    str=sprintf('Time = %5.2f days',t(end));
    titles=sgtitle(str);
    drawnow
else
    % Figure exists, update data
    for i=1:param.Nx
        set(plots(1,i),'XData',t,'YData',x(i,:))
    end
    for i=1:param.Ns
        set(plots(2,i),'XData',t,'YData',S(i,:))
        set(plots(3,i),'XData',t,'YData',bflux(i,:))
        set(plots(4,i),'XData',z,'YData',Sb(i,:))
        %set(plots(4,i),'XData',[z(end),z(end)+param.LL],'YData',[Sb(i,end),S(i,end)]);
    end
    set(plots(5,1),'XData',t,'YData',1E6*thickness)
    str=sprintf('Time = %5.2f days',t(end));
    set(titles,'String',str)
    drawnow
end
end