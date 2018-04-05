% 
% f = openfig('JacobiErrorBig.fig');
% H = findobj(f,'type','line');
% x_dataJ = get(H,'XData');
% y_dataJ = get(H,'YData');
% clf
% 
% 
% 
% f = openfig('GSsmallvbig.fig');
% H = findobj(f,'type','line');
% x_dataGS = get(H,'XData');
% y_dataGS = get(H,'YData');
% clf
% 

% errorplot(length(x_dataJ{1,:}),y_dataJ{3,:},y_dataJ{2,:},y_dataJ{1,:})
% errorplot(length(x_dataGS{7,:}),y_dataGS{7,:},y_dataGS{6,:},y_dataGS{5,:})
% legend({'J Norm 1','J Norm 2', 'J Norm \infty','GS Norm 1','GS Norm 2', 'GS Norm \infty'},'FontSize',20)
% 
% f = openfig('SORBigError.fig');
% H = findobj(f,'type','line');
% x_dataSOR = get(H,'XData');
% y_dataSOR = get(H,'YData');
% clf
% 
% f = openfig('SLORBigError.fig');
% H = findobj(f,'type','line');
% x_dataSLOR = get(H,'XData');
% y_dataSLOR = get(H,'YData');
% clf
% 
% f = openfig('ADIErrorBigFix.fig');
% H = findobj(f,'type','line');
% x_dataADI = get(H,'XData');
% y_dataADI = get(H,'YData');
% clf


% errorplot(length(x_dataSOR{6,:}),y_dataSOR{6,:},y_dataSOR{5,:},y_dataSOR{4,:})
% errorplot(length(x_dataSLOR{6,:}),y_dataSLOR{6,:},y_dataSLOR{5,:},y_dataSLOR{4,:})
% errorplot(length(x_dataADI{1,:}),y_dataADI{3,:},y_dataADI{2,:},y_dataADI{1,:})

legend({'SOR Norm 1','SOR Norm 2', 'SOR Norm \infty','SLOR Norm 1','SLOR Norm 2','SLOR Norm \infty','ADI Norm 1','ADI Norm 2', 'ADI Norm \infty'},'FontSize',20)