function [state]=Props(state)
%Check Dimensions
if state.dimensions == 1
    h=state.enthalpy./state.mass;
else
    h=state.specificenthalpy;
end
%Check for Phase Change Fluid
if strcmp(state.type,'solid')
    M_w=(state.mass.*state.Cp);
    state.temperature=state.enthalpy./M_w;
elseif strcmp(state.type,'2phasefluid')
    p=state.pressure;
    state.temperature=...
        interp2DUG(state.fluid.PP,state.fluid.HH,state.fluid.T_P_H,p,h);
    state.density=...
        interp2DUG(state.fluid.PP,state.fluid.HH,state.fluid.D_P_H,p,h);
    state.quality=...
        interp2DUG(state.fluid.PP,state.fluid.HH,state.fluid.Q_P_H,p,h);
elseif strcmp(state.type,'1phasefluid')
    state.temperature=interp1DUG(state.fluid.HH,state.fluid.T_H,h);
    %state.density=interp1DUG(state.fluid.HH,state.fluid.D_H,h);
end