tic
clear all
% clf
dx = 0.25;
dy = 0.25;
% SOR = 1.2563; %SOR OPT
SOR = 1.24;
x = 0:dx:6;
y = 0:dy:4;
beta = dx/dy;
[X,Y] = meshgrid(x,y);

xlen = length(x);
ylen = length(y);

u = zeros(ylen,xlen);

%BCs
u(1,:) = 1;
u(:,end) = 1;
u(:,1) = 1;
u(end,:) = 1;
u(1,8:18) = 0;
l1norm = 1;
l1norm_target = 1E-6;
iter = 0;
while l1norm > l1norm_target
    for j=2:ylen-1
        utemp = zeros(xlen,xlen);
        b = zeros(xlen,1);
        un = u;
        for i = 2:xlen-1
                utemp(i,i) = 2*(1+beta^2);
                utemp(i,i+1) = -SOR;
                utemp(i,i-1) = -SOR;
                b(i) = 2*(1+beta^2)*(1-SOR)*un(j,i)+SOR*beta^2*(un(j+1,i)+u(j-1,i));
        end
        utemp(1,1) = 1;
        utemp(end,end) = 1;
        b(1) = 1;
        b(end) = 1;
        utemp2 = td(utemp,b);
        u(j,2:end-1) = utemp2(2:end-1);
        
    end
    iter = iter + 1;
%     l1norm(iter) = norm(u-un);
    l1norm(iter) = norm(u(:)-un(:),1);
    l2norm(iter) = norm(u(:)-un(:));
    linfnorm(iter) = norm(u(:)-un(:),Inf);
    if iter>10000
       break 
    end
end

% contour(X,Y,u)
iter
toc
solnplot2(X,Y,u)
errorplot(iter,l1norm,l2norm,linfnorm)