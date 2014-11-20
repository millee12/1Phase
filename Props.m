function [state]=Props(state)
fluid=state.fluid;
%Check Dimensions
if state.dimensions == 1
    h=state.enthalpy./state.mass;
else
    h=state.specificenthalpy;
end
%Check for Phase Change Fluid
if strcmp(state.type,'2phasefluid')
    p=state.pressure;
    state.temperature=interp2DUG(fluid.PP,fluid.HH,fluid.T_P_H,p,h);
    state.density=interp2DUG(fluid.PP,fluid.HH,fluid.D_P_H,p,h);
    state.quality=interp2DUG(fluid.PP,fluid.HH,fluid.Q_P_H,p,h);
elseif strcmp(state.type,'1phasefluid')
    state.temperature=interp1DUG(fluid.HH,fluid.T_H,h);
    state.density=interp1DUG(fluid.HH,fluid.D_H,h);
else
    error('unable to read fluid type in Props')
end