% clear all
% clf
tic
dx = 0.25;
dy = 0.25;
SOR = 1.725; %SOR OPT
% SOR = 1.24; %SOR OPT
x = 0:dx:6;
y = 0:dy:4;
beta = dx/dy;
[xx,yy] = meshgrid(x,y);

xlen = length(x);
ylen = length(y);

u = zeros(ylen,xlen);

%BCs
u(1,:) = 1;
u(:,end) = 1;
u(:,1) = 1;
u(end,:) = 1;
u(1,8:18) = 0;
l1norm_target=1E-6;
l1norm = 1;
iter = 0;
while l1norm > l1norm_target
    for j=2:ylen-1
        un = u;
        for i = 2:xlen-1
            u(j,i) = (1-SOR)*un(j,i)+(SOR/(2*(1+beta^2)))*((beta^2)*(un(j+1,i)+u(j-1,i))+u(j,i+1)+u(j,i-1));
        end
    end
%     error = norm(u-un);
    iter = iter + 1;
%     l1norm(iter) = sum(abs(u(:))-abs(un(:)));
%     l2norm(iter) = (sum((abs(u(:))-abs(un(:))).^2)).^0.5;
%     linfnorm(iter) = max(abs(u(:))-abs(un(:)));
    l1norm(iter) = norm(u(:)-un(:),1);
    l2norm(iter) = norm(u(:)-un(:));
    linfnorm(iter) = norm(u(:)-un(:),Inf);

end
toc
% contour(X,Y,u)
iter
% xx = 0:0.25:6;
% yy = 0:0.25:4;
% [xx,yy] = meshgrid(yy,xx);
% solnplot2(X,Y,u)
% errorplot(iter,l1norm,l2norm,linfnorm);
%%
% clear all
% clf
% tic
dx = 0.025;
dy = 0.025;
% SOR = 1.725; %SOR OPT
% SOR = 1.24; %SOR OPT
SOR = 1.97;
x = 0:dx:6;
y = 0:dy:4;
beta = dx/dy;
[X,Y] = meshgrid(x,y);
u = interp2(xx,yy,u,X,Y);
xlen = length(x);
ylen = length(y);

% u = 0.5*ones(ylen,xlen);

% %BCs
% u(1,:) = 1;
% u(:,end) = 1;
% u(:,1) = 1;
% u(end,:) = 1;
% % u(1,8:18) = 0;
% u(1,62:180) = 0;
for i=1:10
   u(1,61+i) = 1-0.1*i; 
   u(1,180-i) = 1-0.1*i;
end
l1norm_target=1E-6;
l1norm2 = 1;
l2norm2 = 1;
linfnorm2 = 1;
iter = 0;
tic
while l1norm2 > l1norm_target
    for j=2:ylen-1
        un = u;
        for i = 2:xlen-1
            u(j,i) = (1-SOR)*un(j,i)+(SOR/(2*(1+beta^2)))*((beta^2)*(un(j+1,i)+u(j-1,i))+u(j,i+1)+u(j,i-1));
        end
    end
%     error = norm(u-un);
    iter = iter + 1;
%     l1norm(iter) = sum(abs(u(:))-abs(un(:)));
%     l2norm(iter) = (sum((abs(u(:))-abs(un(:))).^2)).^0.5;
%     linfnorm(iter) = max(abs(u(:))-abs(un(:)));
    l1norm2(iter) = norm(u(:)-un(:),1);
    l2norm2(iter) = norm(u(:)-un(:));
    linfnorm2(iter) = norm(u(:)-un(:),Inf);
%     l1norm2(end)
end
toc
% contour(X,Y,u)
iter
% xx = 0:0.25:6;
% yy = 0:0.25:4;
% [xx,yy] = meshgrid(yy,xx);
% solnplot2(X,Y,u)
errorplot(iter,l1norm2,l2norm2,linfnorm2);