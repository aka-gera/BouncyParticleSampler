
function [t_skeleton, x_skeleton, v_skeleton]=ZZ(  x,v,fun,T, refresh_rate)
i1=2;
j1=3;
dim=size(x,1);
I=ones(dim,1);
t = 0;
finished = 0;
x_skeleton = zeros(dim,1);
v_skeleton = zeros(dim,1);
t_skeleton = zeros(1,1); 

rejected_switches = 0;
accepted_switches = 0;
tic; 

Nabla =fun{i1,j1}(x);
NablaNabla=fun{i1,j1+1}(x);
b =  NablaNabla *I; 
a =  v.* Nabla;
Diff_t_switchES_proposed=zeros(dim,1);
for j=1:dim
    Diff_t_switchES_proposed(j) = switchingtime(a(j),b(j));
end
Diff_t_refresh = -log(rand)/refresh_rate;
compter=1;

while finished==0
    argMin= Diff_t_switchES_proposed== min(Diff_t_switchES_proposed);
    Diff_t_switch_proposed=Diff_t_switchES_proposed(argMin);
    Diff_t = min(Diff_t_switch_proposed,Diff_t_refresh);
    
    if compter> T 
        finished = 1;
    end
    x = x + v * Diff_t;
    t = t + Diff_t ;
    a = a + b * Diff_t;
    Nabla = fun{i1,j1}(x);
    switch_rate = v(argMin) * Nabla(argMin);
    if ( Diff_t_switch_proposed < Diff_t_refresh)
        
        proposedSwitchIntensity = a(argMin);
        aTemp0=a(argMin);
        jj=0;
        bTemp= fun{i1,j1+1}(x)*I ;
        while proposedSwitchIntensity < switch_rate
           b=bTemp +I*(-2)^jj;
            aTemp = aTemp0 + b(argMin) * Diff_t;
            proposedSwitchIntensity = aTemp;
            jj=jj+1;
        end
        
        if proposedSwitchIntensity < switch_rate
            disp('ERROR: Switching rate exceeds bound.')
            disp([' simulated rate: ', num2str(proposedSwitchIntensity)])
            disp([' actual switching rate: ', num2str(switch_rate)])
        end
        if rand* proposedSwitchIntensity <= switch_rate
            % reflect
            v(argMin)= -v(argMin); 
            a(argMin) = -switch_rate;
            accepted_switches =accepted_switches+ 1;
        else
            a(argMin) = switch_rate;
            rejected_switches = rejected_switches+1;
        end
        % update refreshment time and switching time bound
        Diff_t_refresh = Diff_t_refresh - Diff_t_switch_proposed;
    else
        % so we refresh;%
        argMin = randi(1:dim);
        a(argMin) = v(argMin) * Nabla(argMin);
        Diff_t_refresh = -log(rand)/refresh_rate;
    end
    for j=1:dim
        Diff_t_switchES_proposed(j) = switchingtime(a(j),b(j));
    end
    x_skeleton(:,compter)=x;
    v_skeleton(:,compter)=v;
    t_skeleton(1,compter)=t;
    compter=compter+1;
    
end
ttime=toc;
nom='ZZ';
disp([nom,'_H: ratio of accepted switches: ', num2str(accepted_switches/(accepted_switches+rejected_switches))])
disp([nom,'_H: number of proposed switches: ', num2str(accepted_switches + rejected_switches)])
disp([nom,'_H: Time of simulation: ', num2str(ttime)])
 
