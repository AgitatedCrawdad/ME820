%MESH SETUP%
nx = 25;
ny = 17;
xmin = 0;
xmax = 6;
ymin = 0;
ymax = 4;

dx = (xmax - xmin) / (nx - 1);
dy = (ymax - ymin) / (ny - 1);

u = zeros(nx,ny);
ud = zeros(nx,ny);

x = linspace(xmin,xmax,nx);
y = linspace(ymin,ymax,ny);

l1norm_target = 1E-6;
l1norm = 1;
%BOUNDARY CONDITIONS%
u(1:7,1) = 1;
u(8:18,1) = 0;
u(19:25,1) = 1;
u(:,end) = 1;
u(1,:) = 1;
u(end,:) = 1;
iter = 0;
tic
while l1norm > l1norm_target
    un = u;
    for i=2:nx-1
    for j=2:ny-1
        u(i,j) = ((un(i+1,j)+un(i-1,j))*dy^2+ (un(i,j+1)+un(i,j-1))*dx^2 )/((dx^2+dy^2)*2);
    end
    end

    iter = iter + 1;
    l1norm(iter) = sum(abs(u(:))-abs(un(:)));
    l2norm(iter) = (sum((abs(u(:))-abs(un(:))).^2)).^0.5;
    linfnorm(iter) = max(abs(u(:))-abs(un(:)));
end
toc
%% PLOTS 
% l2norm(1) = [];
xx = 0:0.25:6;
yy = 0:0.25:4;
xx = 0:dx:6;
yy = 0:dy:4;

[xx,yy] = meshgrid(yy,xx);
contour(yy,xx,u);
