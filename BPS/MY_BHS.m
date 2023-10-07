
function [t_skeleton, x_skeleton, v_skeleton,W]=MY_BHS(  x,v,fun,T, refresh_rate )
Phi=pi/4;
dim=size(x,1); %dimension
iDentity=eye(dim )  ;
t = 0.0;
finished = 0;
x_skeleton = zeros(dim,1);
v_skeleton = zeros(dim,1);
t_skeleton = zeros(1,1);
W=zeros(1,1);
rejected_switches = 0;
accepted_switches = 0;
tic;

i1=2;
ForceGval =fun{i1,1}(x);
NablaForce=fun{i1,2}(x);
aForce =  (v)'  * ForceGval;
bForce = (v)' *NablaForce* v;
Nabla =fun{i1,3}(x);
NablaNabla=fun{i1,4}(x);
aNabla =  (v)'  * Nabla;
bNabla = (v)' *NablaNabla* v;

Diff_t_switch_proposedForce = switchingtime(aForce,bForce);
Diff_t_switch_proposedNabla = switchingtime(aNabla,bNabla);
Diff_t_switch_proposed=min(Diff_t_switch_proposedForce,Diff_t_switch_proposedNabla);
if Diff_t_switch_proposedForce<Diff_t_switch_proposedNabla
    j1=1;
    a=aForce;
    b=bForce;
    Q=NablaForce;
else
    j1=3;
    a=aNabla;
    b=bNabla;
    Q=NablaNabla;
end

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
        finished = 1;
    end
    
    x = x + v * Diff_t/2;
    v=v+(fun{2,1}(x) )*Diff_t;
    x = x + v * Diff_t/2;
    Nabla = fun{i1,j1}(x);
    t = t + Diff_t ;
    a = a + b * Diff_t;
    
    if ( Diff_t_switch_proposed < Diff_t_refresh)
        switch_rate = (v)' * Nabla;
        proposedSwitchIntensity = a;
        
        % Trying to find the best upper bound for the intensity
        jj=0;
        QTemp= fun{i1,j1+1}(x);
        aTemp=a;
        while proposedSwitchIntensity < switch_rate
            Q= QTemp+iDentity*(-1.2)^(jj);
            b=v'*Q*v;
            a  = aTemp + b* Diff_t;
            proposedSwitchIntensity = a;
            jj=jj+1;
        end
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
        
        Diff_t_switch_proposed = switchingtime(a,b);
    else
        % so we refresh;%
        v =cos(Phi)*v+sin(Phi)* randn(dim,1);
        
        ForceGval =fun{i1,1}(x);
        NablaForce=fun{i1,2}(x) ;
        aForce =  (v)'  * ForceGval;
        bForce = (v)' *NablaForce* v ;
        
        Nabla =fun{i1,3}(x);
        aNabla =  (v)'  * Nabla;
        bNabla = (v)' *Q* v;
        
        Diff_t_switch_proposedForce = switchingtime(aForce,bForce);
        Diff_t_switch_proposedNabla = switchingtime(aNabla,bNabla);
        Diff_t_switch_proposed=min(Diff_t_switch_proposedForce,Diff_t_switch_proposedNabla) ;
        if Diff_t_switch_proposedForce<Diff_t_switch_proposedNabla
            j1=1;
            a=aForce;
            b=bForce;
        else
            j1=3;
            a=aNabla;
            b=bNabla;
        end
        %
        % update upcoming event times
        Diff_t_refresh = -log(rand)/refresh_rate;
    end
    
    x_skeleton(:,compter)=x;
    v_skeleton(:,compter)=v;
    t_skeleton(1,compter)=t;
    compter=compter+1;
    
    
end
ttime=toc;
nom='BPS';
disp([nom,'_Global: ratio of accepted switches: ', num2str(accepted_switches/(accepted_switches+rejected_switches))])
disp([nom,'_Global: number of proposed switches: ', num2str(accepted_switches + rejected_switches)])
disp([nom,'_Global: Time of simulation: ', num2str(ttime)])
