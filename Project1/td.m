function u = td( A, c )
%  Solve the  n x n  tridiagonal system for y:
%
%  [ d(1)  a(1)                                  ] [  u(1)  ]   [  c(1)  ]
%  [ b(2)  d(2)  a(2)                            ] [  u(2)  ]   [  c(2)  ]
%  [       b(3)  d(3)  a(3)                      ] [        ]   [        ]
%  [            ...   ...   ...                  ] [  ...   ] = [  ...   ]
%  [                    ...    ...    ...        ] [        ]   [        ]
%  [                        b(n-1) d(n-1) a(n-1) ] [ u(n-1) ]   [ c(n-1) ]
%  [                                 b(n)  d(n)  ] [  u(n)  ]   [  c(n)  ]
%
%  f must be a vector (row or column) of length n
%  a, b, c must be vectors of length n (note that b(1) and c(n) are not used)


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


