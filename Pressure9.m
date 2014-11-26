function state =Pressure9(state)
if strcmp(state.type,'2phasefluid')
    fluid=state.fluid;
    nf=state.N;
    bp=state.geometry.breakpoints;
    ns=length(bp)-1;
    ls=bp(2:ns+1)-bp(1:ns);
    A=state.geometry.crosssections;
    Vs=sum(A.*ls);
    P=state.pressure(1);
    e=Inf;
    m=state.mass;
    h=state.enthalpy./m;
    wmax=50;
    w=1;
    ehist=Inf*ones(wmax,1);
    pspace=fluid.PP(2)-fluid.PP(1);
    tol=1E-6;
    Pind=((P-fluid.PP(1))./pspace)+1;
    PEind=round(Pind)-1;
    PWind=PEind+2;
    %Relaxation Factor
    rf=1;
    while abs(e) > tol && w < wmax
        PE=fluid.PP(PEind);
        PW=fluid.PP(PWind);
        rhoE=interp2DUG(fluid.PP,fluid.HH,fluid.D_P_H,PE*ones(nf,1),h);
        rhoW=interp2DUG(fluid.PP,fluid.HH,fluid.D_P_H,PW*ones(nf,1),h);
        rho=interp2DUG(fluid.PP,fluid.HH,fluid.D_P_H,P*ones(nf,1),h);
        VE=sum(m./rhoE);
        VW=sum(m./rhoW);
        V=sum(m./rho);
        ev=abs(([VE,V,VW]-Vs)/Vs);
        [e,ind]=min(ev);
        ehist(w+1)=e;
        if ind == 2
            dVdP=rf*sum((VW-VE)/pspace);
            while abs(e) > tol && w < wmax
                dV=Vs-V;
                dP=dV/dVdP;
                P=P+dP;
                rho=interp2DUG(fluid.PP,fluid.HH,fluid.D_P_H,P*ones(nf,1),h);
                V=sum(m./rho);
                e=(V-Vs)/Vs;
                ehist(w+1)=e;
                if P < PE
                    PEind=PEind-1;
                    PWind=PWind-1;
                    break
                elseif P > PW
                    PEind=PEind+1;
                    PWind=PWind+1;
                    break
                end
                w=w+1;
            end
        elseif ind==1
            P=PE;
            PEind=PEind-1;
            PWind=PWind-1;
        elseif ind==3
            P=PW;
            PEind=PEind+1;
            PWind=PWind+1;
        end
        w=w+1;
    end
    state.pressure=P*ones(nf,1);
end