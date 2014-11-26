function [state] = firstmesh(state,m0,FLUID,Tamb)
%Assign Fluid Table
R=load('fluids.mat',FLUID);
state.fluid=R.(genvarname(FLUID));
fluid=state.fluid;
%%%%%
nf=state.N;
ls=state.geometry.breakpoints(2:end)-state.geometry.breakpoints(1:end-1);
vt=sum(ls.*state.geometry.crosssections);
vf=vt/nf;
if strcmp(state.type,'2phasefluid')
    state.mass=(m0/nf)*ones(nf,1);
    state.density=(1/vf)*state.mass;
    
elseif strcmp(state.type,'1phasefluid')
    rho0=interp1DUG(fluid.TT,fluid.D_T,Tamb);
    state.density=rho0*ones(nf,1);
    state.mass=state.density.*vf;
else
    disp('error: give intial density or mass of fluid')
end
%%%%%%make mesh%%%%%%%%%%
[state]=remesh(state);
%%%%%%Get Fluid Properties%%%%%%%
if strcmp(state.type,'2phasefluid')
    %(saturation pressure for Tamb)
    state.pressure=ones(nf,1)*...
        interp2(state.fluid.TT,state.fluid.QQ,state.fluid.P_T_Q,Tamb,.5);
    p=state.pressure(1);
    %Saturated Values
    hl=interp2DUG(fluid.PP,fluid.QQ,fluid.U_P_Q,p,0);
    hlatent=interp2DUG(fluid.PP,fluid.QQ,fluid.U_P_Q,p,.99)-hl;
    vl=1/interp2DUG(fluid.PP,fluid.QQ,fluid.D_P_Q,p,0);
    vlatent=(1/interp2DUG(fluid.PP,fluid.QQ,fluid.D_P_Q,p,.99))-vl;
    state.quality=((state.density.^-1)-vl)/(vlatent);
    if max(state.quality) > 1 || min(state.quality) < 0
        disp(max(state.quality))
        error('initial quality value is non-physical')
    end
    %Enthalpy
    state.enthalpy=state.mass.*(hl+state.quality*hlatent);

elseif strcmp(state.type,'1phasefluid')
    state.temperature=Tamb*ones(nf,1);
    %Enthalpy
    state.enthalpy=state.mass.*...
        interp1DUG(fluid.TT,fluid.H_T,state.temperature);
else
    disp('error: give intial density or mass of fluid')
end
    %all other properties
    state=Props(state);
    %initialize values
    state.heatrate=zeros(nf,1);
