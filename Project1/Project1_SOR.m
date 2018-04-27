clear all
clf
nx = 25;
ny = 17;
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

l1norm_target = 0.0001;
l1norm = 1;


u(:,1) = 1;
u(1:7,end) = 0;
u(8:end,end) = 1;
u(1,:) = 1;
u(end,:) = 1;
omega = 1.7;
iter = 0;
while l1norm > l1norm_target
    un = u;
    for i=2:nx-1
    for j=2:ny-1
        u(i,j) = ((un(i+1,j)+u(i-1,j))*dy^2+ (un(i,j+1)+u(i,j-1))*dx^2 )/((dx^2+dy^2)*2);
        u(i,j) = un(i,j) + omega*(u(i,j)-un(i,j));
    end
    end
    u(1:6,1) = 1;
    u(7:19,1) = 0;
    u(20:25,1) = 1;
    u(:,end) = 1;
    u(1,:) = 1;
    u(end,:) = 1;
    l1norm = sum(abs(u(:))-abs(un(:)))/(sum(abs(un(:))));
    iter = iter + 1;
end
xx = 0:0.25:6;
yy = 0:0.25:4;
[xx,yy] = meshgrid(yy,xx);
figure(1)
surface(yy,xx,u);
% colorbar
% shading interp
% figure(2)
% hold on
contour(yy,xx,u)
% grid on
toc
iter
