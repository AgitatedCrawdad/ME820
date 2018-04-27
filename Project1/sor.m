clf
n = 25;
d = -4*ones(n,1); 
a = ones(n,1);
a(end) = 1;
b = ones(n,1);
b(1) = 1;
A = diag(d,0) + diag(b(2:end),-1) + diag(b(1:end-1),1);
c = zeros(n,1);
c(1:5) = -1;

u = zeros(n,1);
% for i = 1:12
%     for j=2:n
%         d(j) = d(j) -(b(j)/(d(j-1)))*a(j-1);
%         c(j) = c(j) -(b(j)/(d(j-1)))*c(j-1);
%     end
%     for k=n-1:-1:1
%        u(k) = (c(k) - a(k)*u(k+1))/d(k);
%     end
% end

% for j=2:n
%     d(j) = d(j) -(b(j)/(d(j-1)))*a(j-1);
%     c(j) = c(j) -(b(j)/(d(j-1)))*c(j-1);
% end
% for k=n-1:-1:1
%    u(k) = (c(k) - a(k)*u(k+1))/d(k);
% end
c1 = c;
[c,u] = td(b,d,a,c);

y = reshape(u,sqrt(n),sqrt(n));
colorbar
surface(y);
shading interp