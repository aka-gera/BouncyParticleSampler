%  Energy

fun = Myfunction1();

T=1000000; 
refresh_rate=1;
x=rand(2,1); % initial position
v=rand(2,1);%  %initial velovity
dim=size(x,1);
%
for cas=1:5
    tic;
    switch cas
        case 1
            Sig= eye(dim);
            mala_stepsize = 0.01/dim^(4/3);
            xMALA=  MALA(x,mala_stepsize,  T,fun);
        case 2
            [~, xBPSL, ~]=BPS_Local(  x,v,fun,T, refresh_rate);
        case 3
            [~, xZZ, ~] =ZZ(  x,v,fun,T, refresh_rate);
        case 4
            [~, xBPSG, ~] =BPS_Global(  x,v,fun,T, refresh_rate );
        case 5
           [~, xBHS, ~] =MY_BHS(  x,v,fun,T, refresh_rate );
    end
    
end
%  plot graph KL and Histogram
% xT=xBPS;
% tic;
LaB={'BPS Global','BPS Local','ZZ','BHS','MALA'}; 
Col={'r+-','b','k','g','r*-'};
xT={xBPSG,xBPSL,xZZ,xBHS,xMALA};
GraphG(fun,xT,LaB,Col)
% toc
% figure(4)   
% clf
% plot([xT(1,1:end-2);xT(1,2:end-1)], [xT(2,1:end-2);xT(2,2:end-1)])




