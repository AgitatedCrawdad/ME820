function solnplot(yy,xx,u)
figure(1)
g=gcf;
g.Units='inches';
g.Position=[-18 1 11.25 7.5];
contour(yy,xx,u,'ShowText','on','LineWidth',3)
% colorbar
set(gca,'fontsize',20)
% colorbar
xlabel('x [m]','FontSize',20)
ylabel('y [m]','FontSize',20)
title('Solution','FontSize',20)
y = linspace(0,4,17);
x = linspace(0,6,25);
for k = 1:length(y)
line([x(1) x(end)], [y(k) y(k)],'Color',[0 0 0],'LineWidth',1);
end
% Vertical grid
for k = 1:length(x)
line([x(k) x(k)], [y(1) y(end)],'Color',[0 0 0],'LineWidth',1);
end