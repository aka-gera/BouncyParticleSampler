 

function [t_skeleton, x_skeleton, v_skeleton,W]=BPSG_H(  x,v,fun, nnn,T, refresh_rate )

dim=size(x,1);
t = 0.0;
finished = 0;
x_skeleton = zeros(dim,1);
v_skeleton = zeros(dim,1);
t_skeleton = zeros(1,1);
W=zeros(1,1);
rejected_switches = 0;
accepted_switches = 0;
tic;

i1=nnn(1);
Nabla =fun{i1,1}(x);
Q=fun{i1,2}(x);
a =  (v)'  * Nabla;
b = (v)' *Q* v;

Diff_t_switch_proposed = switchingtime(a,b);

if (refresh_rate == 0.0)
    Diff_t_refresh = Inf;
else
    Diff_t_refresh = -log(rand)/refresh_rate;
end
compter=1;
ikk=1;
while finished==0
    Diff_t = min(Diff_t_switch_proposed,Diff_t_refresh);
    
    if compter> T
        Diff_t = T - t;
        finished = 1;
    end
    x = x + v * Diff_t;
    Nabla = fun{i1,1}(x);
    t = t + Diff_t ;
    a = a + b * Diff_t;
    
    if ( Diff_t_switch_proposed < Diff_t_refresh)
        switch_rate = (v)' * Nabla;
        proposedSwitchIntensity = a;
        
        % Trying to find the best upper bound for the intensity 
        jj=0;
        QTemp= fun{i1,2}(x)  ;
        aTemp=a;
        bTemp= v' * QTemp*v;
        vNorm=v'*v;
        while proposedSwitchIntensity < switch_rate
            b=bTemp + vNorm*(-2)^jj;
            a  = aTemp + b* Diff_t;
            proposedSwitchIntensity = a;
            jj=jj+1;
        end
        Q=QTemp +diag(ones(dim,1))*(-2)^jj;
         % Done
        
        
        if proposedSwitchIntensity < switch_rate
            disp('ERROR: Switching rate exceeds bound.')
            disp([' simulated rate: ', num2str(proposedSwitchIntensity)])
            disp([' actual switching rate: ', num2str(switch_rate)])
        end
        if rand* proposedSwitchIntensity <= switch_rate
            % reflect
            v = Reflect(Nabla,v);
            a = -switch_rate;
            b = transpose(v) * Q*v;
            accepted_switches =accepted_switches+ 1;
            W(1,ikk)=compter;
            ikk=ikk+1;
        else
            a = switch_rate;
            rejected_switches = rejected_switches+1;
        end
        % update refreshment time and switching time bound
        Diff_t_refresh = Diff_t_refresh - Diff_t_switch_proposed;
    else
        % so we refresh;%
        v = randn(dim,1);
        a = transpose(v) * Nabla;
        b =transpose(v) * Q* v;
        % update upcoming event times
        Diff_t_refresh = -log(rand)/refresh_rate;
    end
    
    x_skeleton(:,compter)=x;
    v_skeleton(:,compter)=v;
    t_skeleton(1,compter)=t;
    compter=compter+1;
    
    Diff_t_switch_proposed = switchingtime(a,b);
    
end
ttime=toc;
nom='BPSG';
disp([nom,'_H: ratio of accepted switches: ', num2str(accepted_switches/(accepted_switches+rejected_switches))])
disp([nom,'_H: number of proposed switches: ', num2str(accepted_switches + rejected_switches)])
disp([nom,'_H: Time of simulation: ', num2str(ttime)])
 
