function plt_structure(POINTS,PANELS,CG,Npoints,Npanels)
% PLOT - data
% this function plots the structure 

figure(1)
for ii = 1:Npanels
    x1 = PANELS(ii).start.coords(1);
    y1 = PANELS(ii).start.coords(2);
    x2 = PANELS(ii).finish.coords(1);
    y2 = PANELS(ii).finish.coords(2);
    quiver(x2,y2,x1-x2,y1-y2,'LineWidth',4);
    hold on
end 

% panels plot 
for ii = 1:Npanels
    x1 = PANELS(ii).start.coords(1);
    y1 = PANELS(ii).start.coords(2);
    x2 = PANELS(ii).finish.coords(1);
    y2 = PANELS(ii).finish.coords(2);
    x = [x1, x2];
    y = [y1, y2];
    plot(x,y,'-r','LineWidth',2.5);
end

% points plot 
for ii = 1:Npoints
    TEXT = append('$\textcircled{', num2str(ii), '}$');
    S = 0.025;
    plot(POINTS(ii).coords(1), POINTS(ii).coords(2),'ok','LineWidth',4);
	text(POINTS(ii).coords(1)+S, POINTS(ii).coords(2)+S, TEXT, 'Interpreter', 'latex', 'FontSize',14);
end

% CG plot
plot(CG(1),CG(2),'ob','LineWidth',8);

% initializing names
quiver_name = strings(Npanels,1);
for ii = 1:Npanels
    Q = append('q', num2str(ii));
    quiver_name(ii) = Q;
end

legend(quiver_name);

xlabel('x');
ylabel('y');
grid on 
grid minor 
axis 'equal'
title('structure');