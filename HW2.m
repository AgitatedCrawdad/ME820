%1st order forward scheme
x = 4;
fpexact = (x*cos(x)-3*sin(x))/(x^4);
% double(fpexact)
% h = 0.0001;
% h = logspace(-5,-1,100);
h = linspace(1E-5,1E-1,100);
f = @(x) sin(x)./(x.^3);

fp1 = (f(x+h)-f(x))./h;
error = abs(fpexact - fp1);
loglog(h,error)

%2nd order central scheme
hold on
fp2 = (f(x+h)-f(x-h))./(2*h);
error2 = abs(fpexact - fp2);
loglog(h,error2,'--');
grid on 
legend('1st Order Forward','2nd Order Central')

xlabel('h')
ylabel('|error|')