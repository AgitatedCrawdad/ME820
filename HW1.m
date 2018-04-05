kh = linspace(0,pi,100);
h=1;
pade = (3*sin(kh)./(h*(2+cos(kh))));
center = (-sin(2*kh)+8*sin(kh))/(6*h);
plot(kh,pade,kh,center,'.',kh,kh,'k--')
legend('4th Order Pade','4th Order Central', 'Exact','Location','northwest')
grid on
title('Modified Wave Number')
xlabel('kh')
ylabel('k^\primeh')
axis([0 pi 0 pi])