rho = readmatrix("toPlot/density.txt");
rho1 = readmatrix("toPlot/densityO2.txt");
rho2 = readmatrix("toPlot/densityO.txt");
pres = readmatrix("toPlot/pressure.txt");
temp = readmatrix("toPlot/temp.txt");
vel_tau = readmatrix("toPlot/velocity_tau.txt");
vel_norm = readmatrix("toPlot/velocity_normal.txt");
vel = readmatrix("toPlot/velocity.txt");

FontSize = 20;
%% velocity plotting
figure ("Position", [0, 0, 2000, 700])
t=tiledlayout(1, 3, "TileSpacing", "compact");
nexttile
hold on
plot(vel_tau(:,1),vel_tau(:,2),'-r','LineWidth', 2.0);
grid on
title("tangent velocity",'FontSize',FontSize);
xlabel('length, m','FontSize',FontSize)
ylabel('velocity, m/s','FontSize',FontSize)
xlim([min(vel_tau(:,1))  max(vel_tau(:,1))])
ylim([min(vel_tau(:,2))  max(vel_tau(:,2))])

nexttile
hold on
plot(vel_norm(:,1), vel_norm(:,2),'-r','LineWidth', 2.0);
grid on
title("normal velocity",'FontSize',FontSize);
xlabel('length, m','FontSize',FontSize)
ylabel('velocity, m/s','FontSize',FontSize)
xlim([min(vel_norm(:,1))  max(vel_norm(:,1))])
ylim([min(vel_norm(:,2))  max(vel_norm(:,2))])

nexttile
hold on
plot(vel(:,1), vel(:,2),'-r','LineWidth', 2.0);
grid on
title("velocity",'FontSize',FontSize);
xlabel('length, m','FontSize',FontSize)
ylabel('velocity, m/s','FontSize',FontSize)
xlim([min(vel(:,1))  max(vel(:,1))])
ylim([min(vel(:,2))  max(vel(:,2))])

saveas(gcf,'fig/all_velocity.fig')
saveas(gcf,'fig/all_velocity.png')

%% density plotting
figure ("Position", [0, 0, 2000, 700])
t=tiledlayout(1, 3, "TileSpacing", "compact");
nexttile
hold on
plot(rho1(:,1),rho1(:,2),'-r','LineWidth', 2.0);
grid on
title("O_2 density",'FontSize',FontSize);
xlabel('length, m','FontSize',FontSize)
ylabel('velocity, kg/m^3','FontSize',FontSize)
xlim([min(rho1(:,1))  max(rho1(:,1))])
ylim([min(rho1(:,2))  max(rho1(:,2))])

nexttile
hold on
plot(rho2(:,1), rho2(:,2),'-r','LineWidth', 2.0);
grid on
title("O density",'FontSize',FontSize);
xlabel('length, m','FontSize',FontSize)
ylabel('velocity, kg/m^3','FontSize',FontSize)
xlim([min(rho2(:,1))  max(rho2(:,1))])
ylim([min(rho2(:,2))  max(rho2(:,2))])

nexttile
hold on
plot(rho(:,1), rho(:,2),'-r','LineWidth', 2.0);
grid on
title("density",'FontSize',FontSize);
xlabel('length, m','FontSize',FontSize)
ylabel('velocity, kg/m^3','FontSize',FontSize)
xlim([min(rho(:,1))  max(rho(:,1))])
ylim([min(rho(:,2))  max(rho(:,2))])

saveas(gcf,'fig/all_density.fig')
saveas(gcf,'fig/all_density.png')

%% temperature plotting
figure ("Position", [0, 0, 1000, 1000])
hold on
plot(temp(:,1),temp(:,2),'-r','LineWidth', 2.0);
grid on
title("temperature",'FontSize',FontSize);
xlabel('length, m','FontSize',FontSize)
ylabel('temperature, K','FontSize',FontSize)
xlim([min(temp(:,1))  max(temp(:,1))])
ylim([min(temp(:,2))  max(temp(:,2))])

saveas(gcf,'fig/temperature.fig')
saveas(gcf,'fig/temperature.png')
%% pressure plotting
figure ("Position", [0, 0, 1000, 1000])
hold on
plot(pres(:,1),pres(:,2),'-r','LineWidth', 2.0);
grid on
title("pressure",'FontSize',FontSize);
xlabel('length, m','FontSize',FontSize)
ylabel('pressure, Pa','FontSize',FontSize)
xlim([min(pres(:,1))  max(pres(:,1))])
% ylim([min(pres(:,2))  max(pres(:,2))])

saveas(gcf,'fig/pressure.fig')
saveas(gcf,'fig/pressure.png')
