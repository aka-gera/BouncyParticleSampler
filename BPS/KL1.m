function ff=KL1(xA,xB,IntA,xa,yb,sa,sb)
% KL divergence
 
N1 = hist2(xA(1,:),xA(2,:),xa,yb);  %; compute the number of the time the solution is in a grid
Fa=N1/sum(sum(N1,1),2)/(sa*sb);  % compute the distribution
xxFa=Fa(IntA);
xxF=xxFa./xB;

ind=1e-8*round(1e8*xxF)~=0; 
ff=sum(xxFa(ind).*log(xxF(ind)));% KL divergence

