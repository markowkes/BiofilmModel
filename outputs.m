function outputs(t,x,S,z,bflux,Sb)
%This function organizes all the outputs into one place so that the desired
%plots can be produced in one easy step

figure(1); clf(1);
subplot(2,2,1)
plot(t,x)
title('Biomass Concentration For Filling/Draining Tank')
ylabel('Grams')
xlabel('Time')
hold on
subplot(2,2,2)
plot(t,S)
title('Substrate Concentrations For Filling/Draining Tank')
ylabel('Grams')
xlabel('Time')
hold on
subplot(2,2,3)
plot(t,bflux)
hold on
subplot(2,2,4)
plot(z,Sb)
title('Substrate Concentration within Biofilm')
ylabel('Grams')
xlabel('Depth of Biofilm [m]')

end