nx = 10; % Number of points in mesh
ny = 10; % Number of points in mesh

u = zeros(nx,ny+1); %Initialize u mesh
un = u;
uf = zeros(nx,ny);

v = zeros(nx+1,ny); %Initialize v mesh
vn = v;
vf = zeros(nx,ny);

p = zeros(nx+1,ny+1); %Initialize p mesh
pn = p;
pf = zeros(nx,ny);

dx = (1/(nx-1));
dy = (1/(ny-1));
beta = dx/dy;
dt = 0.01;

error = 1;
error2 = 1;
Re = 10;

u(:,end) = 1;
u(:,end-1) = 1;

iter = 1;

while error > 1E-4
    
    for i = 2:(nx-1)
        for j = 2:ny
            un(i,j) = u(i,j) + dt*(  (u(i+1,j)*u(i+1,j)-u(i-1,j)*u(i-1,j))/(2*dx)...
                -(((u(i,j)+u(i,j+1))/2)*((v(i,j)+v(i+1,j))/2)-((u(i,j)+u(i,j-1))/2)*((v(i+1,j-1)+v(i,j-1))/2))/dy...
                -(((p(i+1,j)-p(i,j))/dx))...
            +(1/Re)*(((u(i+1,j)-2*u(i,j)+u(i-1,j))/dx^2)+((u(i,j+1)-2*u(i,j)+u(i,j-1))/dy^2))    );  
        end
    end
    
    for j =2:(nx+1)
        un(1,j) = 0;
        un(end,j) = 0;
    end
    
    for i =1:(nx)
        un(i,1) = -un(i,2);
        un(i,end) = 2 - un(i,end-1);
    end
% 
    for i = 2:ny
        for j = 2:(nx-1)
            vn(i,j) = v(i,j) - dt*(  (v(i,j+1)*v(i,j+1)-v(i,j-1)*v(i,j-1))/(2*dy)...
                +(((u(i,j)+u(i,j+1))/2)*((v(i,j)+v(i+1,j))/2)-((u(i-1,j)+u(i-1,j+1))/2)*((v(i,j)+v(i-1,j))/2))*dy...
                +(((p(i,j+1)-p(i,j))/dy))...
            -(1/Re)*(((v(i+1,j)-2*v(i,j)+v(i-1,j))/dx^2)+((v(i,j+1)-2*v(i,j)+v(i,j-1))/dy^2))    );  
        end
    end

    for i =1:(nx+1)
        vn(i,1) = 0;
        vn(i,end) = 0;
    end
    
    for j =2:(nx-1)
        vn(1,j) = -vn(2,j);
        vn(end,j) = - vn(end-1,j);
    end
   
    while error2 > 0.0001
%     for i = 2:nx
%        for j = 2:ny
%           pn(i,j) = (   (-2-2*(beta^2))^-1)*...
%               (-p(i-1,j)-p(i+1,j)-(beta^2)*(p(i,j+1)+p(i,j-1))+...
%               ((dx^2)/dt)*((((un(i,j)+(dt*((p(i+1,j)-p(i,j))/dx)))-(un(i-1,j)+(dt*((p(i,j)-p(i-1,j))/dx))))/dx)+...
%               (((vn(i,j)+(dt*((p(i,j+1)-p(i,j))/dy)))-(vn(i,j-1)+(dt*((p(i,j)-p(i,j-1))/dy))))/dy))  );
%        end 
%     end
    for i = 2:nx
       for j = 2:ny
          pn(i,j) = (   (-2-2*(beta^2))^-1)*...
              ((-p(i-1,j)-p(i+1,j)-(beta^2)*(p(i,j+1)+p(i,j-1)))+...
              ((dx^2)/dt)*((((un(i,j)+(dt*((p(i+1,j)-p(i,j))/dx)))-(un(i-1,j)+(dt*((p(i,j)-p(i-1,j))/dx))))/dx)+...
                           (((vn(i,j)+(dt*((p(i,j+1)-p(i,j))/dy)))-(vn(i,j-1)+(dt*((p(i,j)-p(i,j-1))/dy))))/dy))  );
%           T1 = (-2-2*(beta^2))^-1
%           T2 = (dx^2)/dt
%           T3 = ((((un(i,j)+(dt*((p(i+1,j)-p(i,j))/dx)))-(un(i-1,j)+(dt*((p(i,j)-p(i-1,j))/dx))))/dx))
%           T4 = ((((vn(i,j)+(dt*((p(i,j+1)-p(i,j))/dy)))-(vn(i,j-1)+(dt*((p(i,j)-p(i,j-1))/dy))))/dy))
%           T5 = (-p(i-1,j)-p(i+1,j)-(beta^2)*(p(i,j+1)+p(i,j-1)))
       end 
    end

    
    for i = 1:nx
        pn(i,1) = pn(i,2)-2*vn(i,2)/(Re*dy);
        pn(i,end) = pn(i,end-1)-2*vn(i,end-1)/(Re*dy);
    end  
    
    for j = 1:ny+1
        pn(1,j) = pn(2,j)-2*un(2,j)/(Re*dx);
        pn(end,j) = pn(end-1,j)-2*un(end-1,j)/(Re*dx);
    end  

%     for i = 2:nx
%         pn(i,1) = pn(i,2);
%         pn(i,end) = pn(i,end-1);
%     end
%     
%     for j = 1:(ny+1)
%         pn(1,j) = pn(2,j);
%         pn(end,j) = pn(end-1,j);     
%     end
    
    perror(iter) = norm(p(:)-pn(:),1);
    error2 = perror(end)
    p = pn;

    end
    l1norm(iter) = norm(u(:)-un(:),1)+norm(v(:)-vn(:),1);
    error = l1norm(iter);
    u = un;
    v = vn;
    p = pn;
    iter = iter +1

end

for i = 1:nx
    for j = 1:ny
    uf(i,j) = (u(i,j)+u(i,j+1))/2;
    vf(i,j) = (v(i,j)+v(i+1,j))/2;
    pf(i,j) = (p(i,j)+p(i+1,j)+p(i,j+1)+p(i+1,j+1))/4;
    end
end

xx = 0:dx:1;
yy = 0:dy:1;
[xx,yy] = meshgrid(yy,xx);
% 
figure(1)
colorbar
% surface(yy,xx,uf)
contour(yy,xx,uf);
figure(2)
streamline(stream2(xx,yy,uf,vf,xx,yy))
quiver(xx,yy,uf,vf)
