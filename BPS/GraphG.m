
function GraphG(fun,xT,LaB,Col)


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
 Size2=size(xT,2);  
ft=10; 
figure(1)
clf
set(findall(gcf,'-property','FontSize'),'FontSize',ft)
hold on
for i=1:Size2 
SubPlotF(Size2,2,(i-1)*2+1,xT{i},LaB{i},xa,xb,sa, sb,xbb,xf,yf,ft)
end 


figure(2)
clf 
hold on
for i=1:Size2
[SG,klG]=getKL1(xT{i},FA(IndA),IndA,xa,xb,sa,sb) ; 
loglog( (SG), (klG),Col{i},'linewidth',2 ) 
end
legend(LaB)  
grid on
xlabel('Sample size','Interpreter','tex','FontSize',ft)
ylabel('KL','Interpreter','tex','FontSize',ft)