function = Outputs(t,x,S)
%This function organizes all the outputs into one place so that the desired
%plots can be produced in one easy step

%Call on additional variables from other functions for plots
[Sb,bflux,flux,z]=Diffusion(Lf,So,mumax,Xb,Yxs,De);

figure(3); clf(3);
subplot(3,2,1)
plot(t,x)
title('Biomass Concentration For Filling/Draining Tank')
ylabel('Grams')
xlabel('Time')
hold on
subplot(3,2,2)
plot(t,S)
title('Substrate Concentrations For Filling/Draining Tank')
ylabel('Grams')
xlabel('Time')
hold on
subplot(3,2,3)
plot(t,bflux)
hold on
subplot(3,2,4)
plot(z,Sb)
title('Substrate Concentration within Biofilm')
ylabel('Grams')
xlabel('Depth of Biofilm [m]')

end