clear all
clf
nx = 241;
ny = 161;
xmin = 0;
xmax = 6;
ymin = 0;
ymax = 4;
tic
dx = (xmax - xmin) / (nx - 1);
dy = (ymax - ymin) / (ny - 1);

u = zeros(nx,ny);
ud = zeros(nx,ny);

x = linspace(xmin,xmax,nx);
y = linspace(ymin,ymax,ny);


l1norm_target = 1E-6;
l1norm = 1;

% while l1norm > l1norm_target
%     un = u;
%     for i=2:nx-1
%     for j=2:ny-1
%         u(i,j) = ((un(i+1,j)+un(i-1,j))*dy^2+ (un(i,j+1)+un(i,j-1))*dx^2 )/((dx^2+dy^2)*2);
%     end
%     end
%     u(:,end) = 1;
%     l1norm = sum(abs(u(:))-abs(un(:)))/(sum(abs(un(:))));
% end


% u(:,1) = 1;
% u(1:7,end) = 0;
% u(8:end,end) = 1;
% u(1,:) = 1;
% u(end,:) = 1;
u(1:((length(u(:,1))-1)/4)+1,1)=1;
u(((length(u(:,1))-1)/4)+2:(length(u(:,1))-1)-((length(u(:,1))-1)/4),1)=0;
u((length(u(:,1))-1)-((length(u(:,1))-1)/4)+1:end,1)=1;
for i=1:10
   u(61+i,1) = 1-0.1*i; 
   u(180-i,1) = 1-0.1*i;
end
% % u(1:7,1) = 1;
% % u(8:18,1) = 0;
% % u(19:25,1) = 1;
u(:,end) = 1;
u(1,:) = 1;
u(end,:) = 1;
iter = 0;
uic = u;
while l1norm > l1norm_target
    un = u;
    for i=2:nx-1
    for j=2:ny-1
        u(i,j) = ((un(i+1,j)+un(i-1,j))*dy^2+ (un(i,j+1)+un(i,j-1))*dx^2 )/((dx^2+dy^2)*2);
        
%         u(1:((length(u(:,1))-1)/4)+1,1)=1;
%         u(((length(u(:,1))-1)/4)+2:(length(u(:,1))-1)-((length(u(:,1))-1)/4),1)=0;
%         u((length(u(:,1))-1)-((length(u(:,1))-1)/4)+1:end,1)=1;
%         u(1:7,1) = 1;
%         u(8:18,1) = 0;
%         u(19:25,1) = 1;
%         u(:,end) = 1;
%         u(1,:) = 1;
%         u(end,:) = 1;
    end
    end

    iter = iter + 1;
    l1norm(iter) = sum(abs(u(:))-abs(un(:)));
    l2norm(iter) = (sum((abs(u(:))-abs(un(:))).^2)).^0.5;
    linfnorm(iter) = max(abs(u(:))-abs(un(:)));
%     l1norm = sum(abs(u(:))-abs(un(:)))/(sum(abs(un(:))));
%     l1norm(iter) = norm(u(:)-un(:),1);
%     l2norm(iter+1) = norm(u(:)-un(:));
%     linfnorm(iter) = norm(u(:)-un(:),Inf);
    
end
% PLOTS 
% l2norm(1) = [];
xx = 0:0.25:6;
yy = 0:0.25:4;
xx = 0:dx:6;
yy = 0:dy:4;
toc
iter
[xx,yy] = meshgrid(yy,xx);
solnplot(yy,xx,u);
figure(1)
% surface(yy,xx,u);
g=gcf;
g.Units='inches';
g.Position=[-18 0 11.25 7.5];
contour(yy,xx,u)
% colorbar
xlabel('x [m]')
ylabel('y [m]')
grid on

errorplot(iter,l1norm,l2norm,linfnorm);
figure(2)

iters = linspace(1,iter,iter);
semilogy(iters,l1norm,iters,l2norm,iters,linfnorm)
legend('Norm 1','Norm 2', 'Norm \infty')
xlabel('Number of Iterations')
ylabel('Error')
grid on
g=gcf;
g.Units='inches';
g.Position=[-18 0 11.25 7.5];
% 
% 
figure1=figure(3);

g=gcf;
g.Units='inches';
g.Position=[-18 1 11.25 7.5];
g.Position();
% surface(yy,xx,uic);
% shading interp
axes2 = axes('Parent',figure1); % Create axes
hold(axes2,'on');
imagesc([0 241],[0 161],transpose(uic),'Parent',axes2); % Create image
axes2.YDir='normal';
y = linspace(0,161,161);
x = linspace(0,241,241);
% Horizontal grid
for k = 1:length(y)
line([x(1) x(end)], [y(k) y(k)],'Color',[0 0 0],'LineWidth',1);
end
% Vertical grid
for k = 1:length(x)
line([x(k) x(k)], [y(1) y(end)],'Color',[0 0 0],'LineWidth',1);
end
x1 = 1:20:241;
y1 = 1:20:161;
set(gca,'XTick',x1 );
set(gca,'XTickLabel',0.025*x1-0.025 );
set(gca,'YTick',y1 );
set(gca,'YTickLabel',0.025*y1-0.025 );
xlim([0 241]);
ylim([0 161]);
% title('Mesh','FontSize',20);
xlabel('X Distance [m]','FontSize',20);
ylabel('Y Distance [m]','FontSize',20);
