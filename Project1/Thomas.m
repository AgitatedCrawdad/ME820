function u = td( A, c )

n = length(c);
u = zeros(n,1); 
d = diag(A);
a = diag(A,1);
a(n) = 1;
b(2:n) = diag(A,-1);
b(1) = 1;
for j=2:n
    d(j) = d(j) -(b(j)/(d(j-1)))*a(j-1);
    c(j) = c(j) -(b(j)/(d(j-1)))*c(j-1);
end
u(n) = c(n)/d(n);
for k=n-1:-1:1
   u(k) = (c(k) - a(k)*u(k+1))./d(k);
end


