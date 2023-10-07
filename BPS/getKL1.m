function [SZB,kl1B]=getKL1(xBPS,DF,IndA,xa,xb,sa,sb) 
kl1B=zeros(1);  
zxz=size(xBPS,2);


SZB=1:round(zxz/100):zxz;
for ii=1:length(SZB)
xB=xBPS(:,1:SZB(ii)); 
kl1B(ii)=KL1(xB,DF,IndA,xa,xb,sa,sb) ; 
end
 