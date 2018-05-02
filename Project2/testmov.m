clear
nx = 40; % Number of points in mesh
ny = 40; % Number of points in mesh

u = zeros(nx,ny+1); %Initialize u mesh
un = u;
uf = zeros(nx,ny);
F = u;

v = zeros(nx+1,ny); %Initialize v mesh
vn = v;
vf = zeros(nx,ny);
G = v;

p = zeros(nx+1,ny+1); %Initialize p mesh
pn = p;
pf = zeros(nx,ny);

dx = (1/(nx-1));
dy = (1/(ny-1));
beta = dx/dy;
dt =   0.001;
time = 5.000; %Time in seconds
nloops = time/dt;
% nloops = 250;
error = 1;
error2 = 1;
error2_target = 1E-6;
l1norm_target = 1E-6;
l1norm = 1;
Re = 400;

if dt > Re*0.25*(dx^2)
    disp('dt must be smaller, minimum dt size is...')
    dtsize = Re*0.25*(dx^2)
    return
end

u(2:end-1,end) = 1;
u(2:end-1,end-1) = 1;

iter = 1;

xx = 0:dx:1;
yy = 0:dy:1;
% xxx = 0:2*dx:1;
% yyy = 0:2*dy:1;
[xx,yy] = meshgrid(yy,xx);
% [xxx,yyy] = meshgrid(yyy,xxx);
xxx = 0.5*ones(20,1);
yyy = linspace(0,1,20);
while  iter < nloops
% while l1norm(end) > l1norm_target
    %% Calculates F
    for i = 2:(nx-1)
        for j = 2:ny
            F(i,j) = u(i,j) + dt*(  -((u(i+1,j)*u(i+1,j)-u(i,j)*u(i,j))/(dx))...
                -(((u(i,j)+u(i,j+1))/2)*((v(i,j)+v(i+1,j))/2)-((u(i,j)+u(i,j-1))/2)*((v(i+1,j-1)+v(i,j-1))/2))/dy...
            +(1/Re)*(((u(i+1,j)-2*u(i,j)+u(i-1,j))/dx^2)+((u(i,j+1)-2*u(i,j)+u(i,j-1))/dy^2))    );  
        end
    end
    for j =2:(nx+1)
        F(1,j) = 0;
        F(end,j) = 0;
    end
    
    for i =1:(nx)
%         F(i,1) = 100;
%         F(i,end) = 100;
        F(i,1) = -F(i,2);
        F(i,end) = 2 - F(i,end-1);
    end
    %% Calculates G
    for i = 2:ny
        for j = 2:(nx-1)
            G(i,j) = v(i,j) + dt*(  -((v(i,j+1)*v(i,j+1)-v(i,j)*v(i,j))/(dy))...
                -(((u(i,j)+u(i,j+1))/2)*((v(i,j)+v(i+1,j))/2)-((u(i-1,j)+u(i-1,j+1))/2)*((v(i,j)+v(i-1,j))/2))*dx...
            +(1/Re)*(((v(i+1,j)-2*v(i,j)+v(i-1,j))/dx^2)+((v(i,j+1)-2*v(i,j)+v(i,j-1))/dy^2))    );  
        end
    end
    
    for i =1:(nx+1)
        G(i,1) = 0;
        G(i,end) = 0;
    end
    
    for j =2:(nx-1)
        G(1,j) = -G(2,j);
        G(end,j) = - G(end-1,j);
    end
    
    %% Pressure Poisson Equation Solver
    iter2 = 0;
    while error2 > error2_target

        for i = 2:nx
           for j = 2:ny
              pn(i,j) = (   (-2-2*(beta^2))^-1)*...
                  ((-p(i-1,j)-p(i+1,j)-p(i,j+1)-p(i,j-1))+...
                  ((dx^2)/dt)*((F(i,j)-F(i-1,j))/dx)+((G(i,j)-G(i,j-1))/dy));
           end 
        end
        
%         for i = 1:nx+1
%             pn(i,1) = pn(i,2);
% %             pn(i,end) = 0;
% %             pn(i,end) = pn(i,end-1);
%             pn(i,end) = pn(i,end-1)-2*vn(i,end)/(Re*dy);
%         end  
% 
%         for j = 1:ny+1
%             pn(1,j) = pn(2,j);
%             pn(end,j) = pn(end-1,j);
%         end  


        % BCs
        for i = 1:nx+1
%             pn(i,1) = pn(i,2)-2*vn(i,1)/(Re*dy);
            pn(i,1) = pn(i,2);
%             pn(i,end) = 0;
            pn(i,end) = pn(i,end-1);
%             pn(i,end) = pn(i,end-1)-2*vn(i,end)/(Re*dy);
        end  

        for j = 1:ny+1
            pn(1,j) = pn(2,j);
            pn(end,j) = pn(end-1,j);
%             pn(1,j) = pn(2,j)-2*un(1,j)/(Re*dx);
%             pn(end,j) = pn(end-1,j)-2*un(end,j)/(Re*dx);
        end  

        % Error
        perror(iter2+1) = norm(p(:)-pn(:),1);
        error2 = perror(iter2+1);

        p = pn;
        iter2 = iter2+1;

    end
    
    
    %% Calculates u velocity
    for i = 2:(nx-1)
        for j = 2:ny
            un(i,j) = F(i,j)-dt*(((pn(i+1,j)-pn(i,j))/dx));
        end
    end
    
    for i =1:(nx)
        un(i,1) = -un(i,2);
        un(i,end) = 2 - un(i,end-1);
    end
    for j =2:(nx+1)
        un(1,j) = 0;
        un(end,j) = 0;
    end
    %% Calculates v velocity
    
    for i = 2:ny
        for j = 2:(nx-1)
            vn(i,j) = G(i,j) -dt*(((pn(i,j+1)-pn(i,j))/dy));

        end
    end

 
    for j =2:(nx-1)
        vn(1,j) = -vn(2,j);
        vn(end,j) = - vn(end-1,j);
    end
    for i =1:(nx+1)
        vn(i,1) = 0;
        vn(i,end) = 0;
    end
    
    
    
    l1norm(iter) = norm(u(:)-un(:),1)+norm(v(:)-vn(:),1);
    error = l1norm(iter);
    u = un;
    v = vn;
    p = pn;

    error2 = 1;
    l1norm(end)
%     output = num2str(iter)

    for i = 1:nx
        for j = 1:ny
        uf(i,j,iter) = (u(i,j)+u(i,j+1))/2;
        vf(i,j,iter) = (v(i,j)+v(i+1,j))/2;
        pf(i,j,iter) = (p(i,j)+p(i+1,j)+p(i,j+1)+p(i+1,j+1))/4;
        end
    end
    
    figure(1)
    streamline(stream2(xx,yy,uf(:,:,iter)',vf(:,:,iter)',xx(2:4:end-1,2:4:end-1),yy(2:4:end-1,2:4:end-1)))
    title(strcat('t = ', num2str(dt*iter)));
    xlabel('X-Distance');
    ylabel('Y-Distance');
    MM(iter) = getframe(gcf);
    iter = iter +1;
    clf
end

% figure(1)
% streamline(stream2(xx,yy,uf',vf',xx(2:4:end-1,2:4:end-1),yy(2:4:end-1,2:4:end-1)))
% title(strcat('t = ', num2str(dt*iter)));
% xlabel('X-Distance');
% ylabel('Y-Distance');
% MM = getframe(gcf);

for i = 1:25:5000
    streamline(stream2(xx,yy,uf(:,:,i)',vf(:,:,i)',xx(2:4:end-1,2:4:end-1),yy(2:4:end-1,2:4:end-1)))
%     title(strcat('t = ', num2str(dt*i)));
%     xlabel('X-Distance');
%     ylabel('Y-Distance');
    MM(((i-1)/25)+1) = getframe(gcf);
    clf
end
% MMM = repelem(MM,5);
% movie(MM(2:10).cdata);
% myvideo = VideoWriter('Re100_a.avi');
% open(myvideo)
% writeVideo(myvideo, MM(2:end));
% close(myvideo)
% 
% for i = 1:nx
%     for j = 1:ny
%     uf(i,j) = (u(i,j)+u(i,j+1))/2;
%     vf(i,j) = (v(i,j)+v(i+1,j))/2;
%     pf(i,j) = (p(i,j)+p(i+1,j)+p(i,j+1)+p(i+1,j+1))/4;
%     end
% end

% xx = 0:dx:1;
% yy = 0:dy:1;
% % xxx = 0:2*dx:1;
% % yyy = 0:2*dy:1;
% [xx,yy] = meshgrid(yy,xx);
% % [xxx,yyy] = meshgrid(yyy,xxx);
% xxx = 0.5*ones(20,1);
% yyy = linspace(0,1,20);
% % 
% figure(1)
% colorbar
% surface(yy,xx,pf)
% % contour(yy,xx,pf);
% shading interp
% figure(2)
% streamline(stream2(xx,yy,uf',vf',xxx,yyy))
% hold on
% streamline(stream2(xx,yy,uf',vf',yyy,xxx))
% % % quiver(xx,yy,uf,vf)
% figure(3)
% quiver(yy,xx,uf,vf,0.6,'k-')
% axis equal
% axis([0 1 0 1])
% hold on 
% figure(4)
% pcolor(yy,xx,pf);
% colormap(jet)
% colorbar
% shading interp
% axis equal
% axis([0 1 0 1])
