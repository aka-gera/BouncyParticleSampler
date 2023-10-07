function SubPlotF(mi,ni,k,xM,nom,xa,xb,sa, sb,xbb,xf,yf,ft)
x1=(xM(1,:)); x2=(xM(2,:));
N1 = histcounts(x1,xa);
N2 = histcounts(x2,xb); 
D1=N1/sum(N1)/sa; 
D2=N2/sum(N2)/sb;  

subplot(mi,ni,k)
hold on

bar(xbb,D1)
plot(xa,xf,'r','linewidth',2 )
axis('tight')
grid on
xlabel('x','Interpreter','tex','FontSize',ft)
ylabel('Dist','Interpreter','tex','FontSize',ft)
title(['x vs Dist of  ',nom],'Interpreter','tex','FontSize',ft) 
hold off


subplot(mi,ni,k+1)
hold on

bar( xbb,D2)
plot(xb,yf,'r','linewidth',2 )
axis('tight')
grid on
xlabel('y','Interpreter','tex','FontSize',ft)
ylabel('Dist','Interpreter','tex','FontSize',ft)
title(['y vs Dist of ',nom],'Interpreter','tex','FontSize',ft) 
hold off

