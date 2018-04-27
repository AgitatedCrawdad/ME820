% % 
% f = openfig('JBIGIC.fig');
% H = findobj(f,'type','line');
% x_dataJ = get(H,'XData');
% y_dataJ = get(H,'YData');
% close(f)
% % 
% % 
% % 
% f = openfig('GSBIGIC.fig');
% H = findobj(f,'type','line');
% x_dataGS = get(H,'XData');
% y_dataGS = get(H,'YData');
% close(f)
% % 

% errorplot(length(x_dataJ{1,:}),y_dataJ{3,:},y_dataJ{2,:},y_dataJ{1,:})
% errorplot(length(x_dataGS{1,:}),y_dataGS{3,:},y_dataGS{2,:},y_dataGS{1,:})
% legend({'J Norm 1','J Norm 2', 'J Norm \infty','GS Norm 1','GS Norm 2', 'GS Norm \infty'},'FontSize',20)
% % 
% f = openfig('SORBIGIC.fig');
% H = findobj(f,'type','line');
% x_dataSOR = get(H,'XData');
% y_dataSOR = get(H,'YData');
% close(f)
% % 
% f = openfig('SLORBIGIC.fig');
% H = findobj(f,'type','line');
% x_dataSLOR = get(H,'XData');
% y_dataSLOR = get(H,'YData');
% close(f)
% % 
% f = openfig('ADIBIGIC.fig');
% H = findobj(f,'type','line');
% x_dataADI = get(H,'XData');
% y_dataADI = get(H,'YData');
% close(f)


% errorplot(length(x_dataSOR{1,:}),y_dataSOR{3,:},y_dataSOR{2,:},y_dataSOR{1,:})
% errorplot(length(x_dataSLOR{1,:}),y_dataSLOR{3,:},y_dataSLOR{2,:},y_dataSLOR{1,:})
% errorplot(length(x_dataADI{1,:}),y_dataADI{3,:},y_dataADI{2,:},y_dataADI{1,:})

legend({'SOR Norm 1','SOR Norm 2', 'SOR Norm \infty','SLOR Norm 1','SLOR Norm 2','SLOR Norm \infty','ADI Norm 1','ADI Norm 2', 'ADI Norm \infty'},'FontSize',20)