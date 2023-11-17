vToro = readmatrix("D:\study\course_work\plot\SodaPlotOnly\toPlot\Toro\Velocity.csv");
vHLLE = readmatrix("D:\study\course_work\plot\SodaPlotOnly\toPlot\HLLE\velocity.txt");
vExac = readmatrix("D:\study\course_work\plot\SodaPlotOnly\toPlot\Exac\velocity.txt");
vHLLC = readmatrix("D:\study\course_work\plot\SodaPlotOnly\toPlot\HLLC\velocity.txt");
vToro_new = interp1(vToro(:,1),vToro(:,2),linspace(vToro(1,1),vToro(end,1),100));
clf; hold on
legend;grid on;
plot(1:100,vToro_new,'DisplayName', 'Toro','LineWidth',1.3,'Color','red','LineStyle','none','Marker','o','MarkerSize',3)
plot(linspace(1,100,100),vHLLE(:,2),'DisplayName', 'HLLE','Color','blue','LineWidth',2)
plot(linspace(1,100,200),vExac(:,2),'DisplayName', 'Exac. Riemann','LineWidth',2,'Color','green','LineStyle','-')
plot(vHLLC(:,2),'DisplayName', 'HLLC','LineWidth',2,'Color','magenta','LineStyle','--')
xlabel("y")
ylabel("velocity, m/s")
hold off
saveas(gcf,'velocity.fig')