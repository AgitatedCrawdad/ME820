function errorplot(iter,l1norm,l2norm,linfnorm)
figure(2)
iters = linspace(1,iter,iter);
color = 'c';
semilogy(iters,l1norm,horzcat(color,'-'),'LineWidth',3)
hold on
semilogy(iters,l2norm,horzcat(color,'--'),'LineWidth',3)
semilogy(iters,linfnorm,horzcat(color,':'),'LineWidth',3)
set(gca,'fontsize',12)
legend({'Norm 1','Norm 2', 'Norm \infty'},'FontSize',20)
xlabel('Number of Iterations','FontSize',20)
ylabel('Error','FontSize',20)
title('Error versus Iterations','FontSize',20)
grid on
g=gcf;
g.Units='inches';
g.Position=[-18 1 11.25 7.5];
