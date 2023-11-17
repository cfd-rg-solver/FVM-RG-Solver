vToro = readmatrix("D:\study\course_work\plot\SodaPlotOnly\toPlot\Toro\Velocity.csv");
vHLLE = readmatrix("D:\study\course_work\plot\SodaPlotOnly\toPlot\HLLE\velocity.txt");
vExac = readmatrix("D:\study\course_work\plot\SodaPlotOnly\toPlot\Exac\velocity.txt");
vHLLC = readmatrix("D:\study\course_work\plot\SodaPlotOnly\toPlot\HLLC\velocity.txt");
vToro_new = interp1(vToro(:,1),vToro(:,2),linspace(vToro(1,1),vToro(end,1),100));
clf; hold on
legend;
plot(1:100,vToro_new,'DisplayName', 'Toro','Color','red','LineWidth',1,'LineStyle','-')
plot(vHLLE(:,2),'DisplayName', 'HLLE','Color','blue')
plot(vExac(:,2),'DisplayName', 'Exac. Riemann','Color','green','LineWidth',1.5,'LineStyle','-')
plot(vHLLC(:,2),'DisplayName', 'HLLC','Color','magenta','LineWidth',1.5,'LineStyle','--')
xlabel("y")
ylabel("velocity, m/s")
hold off
saveas(gcf,'velocity.fig')