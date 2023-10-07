%  Energy

fun = Myfunction1();

T=5000; 
refresh_rate=1;
dim=2;
x=rand(dim,1); % initial position
v=rand(dim,1);%  %initial velovity 
%
for cas=5
    tic;
    switch cas
        case 1
            Sig= eye(dim);
            mala_stepsize = 0.01/dim^(4/3);
            xBPS=  MALA(x,mala_stepsize,  T,fun);
        case 2
            [~, xBPS, ~]=BPS_Local(  x,v,fun,T, refresh_rate);
        case 3
            [~, xBPS, ~] =ZZ(  x,v,fun,T, refresh_rate);
        case 4
            [~, xBPS, ~] =BPS_Global(  x,v,fun,T, refresh_rate );
        case 5
           [~, xBPS, ~] =MY_BHS(  x,v,fun,T, refresh_rate );
    end
    
end
%  plot graph KL and Histogram
xT=xBPS;
% tic;
Graph1(fun,xT)
% toc
% figure(4)   
% clf
% plot([xT(1,1:end-2);xT(1,2:end-1)], [xT(2,1:end-2);xT(2,2:end-1)])




