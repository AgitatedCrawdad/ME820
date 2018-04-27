nx = 10; % Number of points in mesh
ny = 10; % Number of points in mesh

u = zeros(nx,ny+1); %Initialize u mesh
un = u;
uf = zeros(nx,ny);

v = zeros(nx+1,ny); %Initialize v mesh
vn = v;
vf = zeros(nx,ny);

p = ones(nx+1,ny+1); %Initialize p mesh
pn = p;
pf = zeros(nx,ny);

dx = (1/(nx-1));
dy = (1/(ny-1));
dt = 0.001;

del = 1.0;
error = 1;
Re = 100;

u(:,end) = 1;
u(:,end-1) = 1;

iter = 1;

while error > 1E-4
    
    for i = 2:(nx-1)
        for j = 2:ny
            un(i,j) = u(i,j) - dt*(  (u(i+1,j)*u(i+1,j)-u(i-1,j)*u(i-1,j))/(2*dx)...
                +(((u(i,j)+u(i,j+1))/2)*((v(i,j)+v(i+1,j))/2)-((u(i,j)+u(i,j-1))/2)*((v(i+1,j-1)+v(i,j-1))/2))/dy...
                +(((p(i+1,j)-p(i,j))/dx))...
            -(1/Re)*(((u(i+1,j)-2*u(i,j)+u(i-1,j))/dx^2)+((u(i,j+1)-2*u(i,j)+u(i,j-1))/dy^2))    );  
        end
    end
    
    for j =2:(nx)
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
    
    for j =1:(nx-1)
        vn(1,j) = -vn(2,j);
        vn(end,j) = - vn(end-1,j);
    end
    
    for i = 2:nx
       for j = 2:ny
          pn(i,j) = p(i,j)-dt*del*(  ( un(i,j)-un(i-1,j) )/dx + (vn(i,j)-vn(i,j-1))/dy );
       end 
    end
    
    for i = 2:nx
        pn(i,1) = pn(i,2);
        pn(i,end) = pn(i,end-1);
    end
    
    for j = 1:(ny+1)
        pn(1,j) = pn(1,j);
        pn(end,j) = pn(end-1,j);     
    end

    l1norm(iter) = norm(u(:)-un(:),1)+norm(v(:)-vn(:),1);
    error(iter) = l1norm(iter);
    u = un;
    v = vn;
    p = pn;
    iter = iter +1;
    error(end)
end

for i = 1:nx
    for j = 1:ny
    uf(i,j) = (u(i,j)+u(i,j+1))/2;
    vf(i,j) = (v(i,j)+v(i+1,j))/2;
    pf(i,j) = (p(i,j)+p(i+1,j)+p(i,j+1)+p(i+1,j+1))/2;
    end
end

xx = 0:dx:1;
yy = 0:dy:1;
[xx,yy] = meshgrid(yy,xx);
% 
colorbar
surface(yy,xx,uf);
% contour(yy,xx,uf);
% streamline(stream2(xx,yy,uf,vf,xx,yy))
