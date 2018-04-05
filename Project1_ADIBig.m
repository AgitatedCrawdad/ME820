clear all
% clf
tic
dx = 0.025;
dy = 0.025;
SOR = 1.295; %SOR OPT
% SOR = 1.24; 
x = 0:dx:6;
y = 0:dy:4;
beta = dx/dy;
[X,Y] = meshgrid(x,y);

xlen = length(x);
ylen = length(y);

u = 0.5*ones(ylen,xlen);

%BCs
u(1,:) = 1;
u(:,end) = 1;
u(:,1) = 1;
u(end,:) = 1;
u(1,62:180) = 0;
for i=1:10
   u(1,61+i) = 1-0.1*i; 
   u(1,180-i) = 1-0.1*i;
end
l1norm_target = 1E-6;
l1norm = 1;
iter = 0;
while l1norm > l1norm_target
%     for j=2:ylen-1
%         utemp = zeros(xlen,xlen);
%         b = zeros(xlen,1);
%         un = u;
%         for i = 2:xlen-1
%                 utemp(i,i) = 2*(1+beta^2);
%                 utemp(i,i+1) = -SOR;
%                 utemp(i,i-1) = -SOR;
%                 b(i) = 2*(1+beta^2)*(1-SOR)*un(j,i)+SOR*beta^2*(un(j+1,i)+u(j-1,i));
%         end
%         utemp(1,1) = 1;
%         utemp(end,end) = 1;
%         b(1) = 1;
%         b(end) = 1;
%         utemp2 = td(utemp,b);
%         u(j,2:end-1) = utemp2(2:end-1);
%         
%     end
    
    for i=2:xlen-1
        utemp = zeros(ylen,ylen);
        b = zeros(ylen,1);
        un = u;
        for j = 2:ylen-1
                utemp(j,j) = 2*(1+beta^2);
                utemp(j,j+1) = -SOR*beta^2;
                utemp(j,j-1) = -SOR*beta^2;
                b(j) = 2*(1+beta^2)*(1-SOR)*un(j,i)+SOR*beta^2*(u(j,i+1)+u(j,i-1));
        end
        
        if(i>61 && i<180)
            utemp(1,1) = 1;
            utemp(end,end) = 1;
            b(1) = 0;
            if i<72
                b(1) = 1-0.1*(72-i);
            end
            if i>169
                b(1) = 0.1*(i-169);
            end
            b(end) = 1;
        else
            utemp(1,1) = 1;
            utemp(end,end) = 1;
            b(1) = 1;
            b(end) = 1;        
        end


        utemp2 = td(utemp,b);
        u(2:end-1,i) = utemp2(2:end-1);
    
    end
    iter = iter + 1;
%     error(iter) = norm(u-un);
    l1norm(iter) = norm(u(:)-un(:),1);
    l2norm(iter) = norm(u(:)-un(:));
    linfnorm(iter) = norm(u(:)-un(:),Inf);
    l1norm(end)
end
toc
% contour(X,Y,u)
iter
solnplot2(X,Y,u);
errorplot(iter,l1norm,l2norm,linfnorm);
