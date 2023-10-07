
function Graph1(fun,xBPS)

%  y = fun{2}(2) % evaluates the second function
am=5;
bm=5;
sa=.1;
sb=.1;
xa=-am:sa:am;
xb=-bm:sb:bm;
xbb=xb(1:end-1);




funX = @(x) exp(- fun{1,1}(x));
Ix = integral(funX,-Inf,Inf) ;% normalization on X dist
funY = @(y) exp(- fun{1,2}(y) );
Iy = integral(funY,-Inf,Inf) ;  %normalization on Y dist

xf=funX(xa) /Ix; %Distribution X
yf=funY(xb)/Iy;% Distribution Y

[xx,yy]=meshgrid(xa,xb);
FA=funX(xx).*funY(yy)/(Ix*Iy); %Distribution XY

IndA=1e-8*round(1e8*FA)~=0;
[SZB,kl1B]=getKL1(xBPS,FA(IndA),IndA,xa,xb,sa,sb) ;
 
mi=1;
ni=2;
ft=10;
figure(1)
clf
set(findall(gcf,'-property','FontSize'),'FontSize',ft)
SubPlotF(mi,ni,1,xBPS,'BPS',xa,xb,sa, sb,xbb,xf,yf,ft)


figure(2)
clf
loglog(SZB,kl1B,'r','linewidth',2 )
axis('tight')
grid on
xlabel('Sample size','Interpreter','tex','FontSize',ft)
ylabel('KL','Interpreter','tex','FontSize',ft)