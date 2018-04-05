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
v1 = ones((nx-2)*(ny-2)-1,1);
v4 = -4*ones((nx-2)*(ny-2),1);
A1 = diag(v1,1);
A2 = diag(v1,-1);
A3 = diag(v4);
A = A1+A2+A3;
b = zeros((nx-2)*(ny-2),1);
b(1:23) = -1;
% b(322:345) = -1;
% for i = 1:13;
%     b(i*23+1) = -1;
% end
iter = 0;
while l1norm > l1norm_target || iter>1000
    un = u;
    w = A\b;
%     for i=2:nx-1
%     for j=2:ny-1
%         u(i,j) = ((un(i+1,j)+u(i-1,j))*dy^2+ (un(i,j+1)+u(i,j-1))*dx^2 )/((dx^2+dy^2)*2);
%     end
%     end
    l1norm = sum(abs(u(:))-abs(un(:)))/(sum(abs(un(:))));
    iter = iter + 1;
end
xx = 0:0.25:6;
yy = 0:0.25:4;
u = reshape(w,23,15);
temp = zeros(25,17);
temp(2:24, 2:16) = u;
% temp(:,1) = 1;
% temp(:,end) = 1;
[xx,yy] = meshgrid(yy,xx);
% figure(1)
surface(yy,xx,temp);
% colorbar
% shading interp
% figure(2)
% hold on
% contour(yy,xx,u)
% grid on
toc
iter
