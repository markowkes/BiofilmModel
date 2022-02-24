function plotSolution(t,X,S,Pb,Sb,Lf,param)
% Outputs - 100 plots with time
    for n=length(t)

        % Update grid
        grid.z  = linspace(0,Lf(n),param.Nz);
        grid.dz = grid.z(2) - grid.z(1);

        figure(1); clf(1)
        subplot(2,3,1)
        plot(t(1:n),X(1:n))
        xlabel('Time')
        ylabel('Tank Particulate Concentrations')
        set(gca,'Fontsize',16)

        subplot(2,3,2)
        plot(t(1:n),S(1:n))
        xlabel('Time')
        ylabel('Tank Substrate Concentrations')
        set(gca,'Fontsize',16)

        subplot(2,3,3)
        plot(t(1:n),Lf(1:n)*1e6)
        xlabel('Time')
        ylabel('Biofilm Thickness (um)')
        set(gca,'Fontsize',16)

        subplot(2,3,4)
        plot(grid.z,Pb(n,:))
        xlabel('Location in Biofilm')
        ylabel('Biofilm Particulate Vol. Fraction')
        set(gca,'Fontsize',16)

        subplot(2,3,5)
        plot(grid.z,Sb(n,:))
        xlabel('Location in Biofilm')
        ylabel('Biofilm Substrate Concentrations')
        set(gca,'Fontsize',16)

    end