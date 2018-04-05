for z = 1:17
    % clear all
    % clf
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
    b = zeros(nx,ny);
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


    u(:,1) = 1;
    u(1:7,end) = 0;
    u(8:end,end) = 1;
    u(1,:) = 1;
    u(end,:) = 1;
    u(3,5:13) = 0;
    % b(3,5:10) = 0.3;
    % b(4:7,14) = 1;
    iter = 0;
    value = 0.0;

    % b(13,4:10) = value;
    b(13,z) = value;
    % b(8:18,end-1) = value;
    % b(3,5:13) = value;
    % for w=1:4
    %     b(w+3,w+9) = value;
    %     b(w+3,9-w) = value;
    % end
    % b(10:15,13) = value;
    % b(10:15,9) = value;
    % b(10:15,5) = value;
    % b(10,9:13) = value;
    % b(15,5:9) = value;
    % 
    % b(18,5:13) = value;
    % b(18:22,5) = value;
    % b(22,5:13) = value;


    while l1norm > l1norm_target
        un = u;
        for i=2:nx-1
        for j=2:ny-1
            u(i,j) = ((un(i+1,j)+u(i-1,j))*dy^2+ (un(i,j+1)+u(i,j-1))*dx^2 -b(i,j)*(dx^2+dy^2))/((dx^2+dy^2)*2);
        end
        end

    %     u(3,5:13) = value;
    %     for w=1:4
    %         u(w+3,w+9) = value;
    %         u(w+3,9-w) = value;
    %     end
    %     u(10:15,13) = value;
    %     u(10:15,9) = value;
    %     u(10:15,5) = value;
    %     u(10,9:13)=value;
    %     u(15,5:9)=value;
    %     
    %     u(18,5:13) = value;
    %     u(18:22,5) = value;
    %     u(22,5:13) = value;
    % %     
    %     u(1:8,1) = 1;
    %     u(8:18,1) = 1;
    %     u(19:25,1) = 1;
    %     u(:,end) = 1;
    %     u(1,:) = 1;
    %     u(end,:) = 1;
    %     

        u(1:8,1) = 1;
        u(1:12,z) =1;
        u(14:end,z) =1;
        
%         u(8:18,1) = 0;
        u(8:18,1) = 0;
        
        
%         u(13,1:12) = 0;
        u(19:25,1) = 1;
        u(:,end) = 1;
        u(8:18,end) = 0;
        u(1,:) = 1;
        u(end,:) = 1;
        iter = iter + 1;
    %     l1norm = sum(abs(u(:))-abs(un(:)))/(sum(abs(un(:))));
    %     l1norm(iter) = sum(abs(u(:))-abs(un(:)));
    %     l2norm(iter) = (sum((abs(u(:))-abs(un(:))).^2)).^0.5;
    %     linfnorm(iter) = max(abs(u(:))-abs(un(:)));
        l1norm(iter) = norm(u(:)-un(:),1);
        l2norm(iter) = norm(u(:)-un(:));
        linfnorm(iter) = norm(u(:)-un(:),Inf);
    end
    xx = 0:0.25:6;
    yy = 0:0.25:4;
    [xx,yy] = meshgrid(yy,xx);
    % figure(1)
    % surface(yy,xx,u);
    % colorbar
    % shading interp
    % g=gcf;
    % g.Units='inches';
    % g.Position=[-18 1 11.25 7.5];

    % figure(2)
    % hold on
    % contour(yy,xx,u)
    % grid on
    toc
    iter
    solnplot(yy,xx,u);
    % errorplot(iter,l1norm,l2norm,linfnorm)
    MM(z) = getframe(gcf);
end
MMM = repelem(MM,20);
movie(MMM);
myvideo = VideoWriter('slot.avi');
open(myvideo)
writeVideo(myvideo, MMM);
close(myvideo)