function [U]=U_calc3(state,regime,mdot_comp,mdot_eev)
persistent mdot_0 t_0
global t
tau=90;
[n,~]=size(state.mass);
U=zeros(n,1);
switch regime
    case 'boiling'
        if mdot_comp == 0
            if mdot_eev==0
                if isempty(t_0)
                   t_0=t; 
                end
                U(:)=0.25/3600+exp((t_0-t)/tau)*(2.5/3600-0.25/3600);
                U(state.quality>=0.99)=...
                    0.01/3600+exp((t_0-t)/tau)*(0.048/3600-0.01/3600);
                U(state.quality<=0.05)=...
                    0.1/3600+exp((t_0-t)/tau)*(0.25/3600-0.1/3600);
            else
                if isempty(mdot_0)
                    mdot_0=mdot_eev;
                end
                U(:)=0.25/3600+(mdot_eev/mdot_0)*(2.5/3600-0.25/3600);
                U(state.quality>=.99)=...
                    0.01/3600+(mdot_eev/mdot_0)*(0.048/3600-0.01/3600);
                U(state.quality<=0.05)=...
                    0.01/3600+(mdot_eev/mdot_0)*(0.25/3600-0.01/3600);
            end
        else
            U(:)=2.5/3600;
            U(state.quality>=.99)=0.048/3600;
            U(state.quality<=0.05)=0.25/3600;
        clear('t_0','mdot_0')
        end
    case 'condensing'
        U(:)=1.2/3600;
        U(state.quality>=.99)=0.388/3600;
        U(state.quality<=0.05)=0.5/3600;
    case 'E-Air'
        U(:)=0.0696/3600;
    case 'C-Air'
        U(:)=0.154/3600;
    case 'none'
        U(:)=0;
    otherwise
        error('Specify Regime')
end