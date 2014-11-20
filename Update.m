function [state,outlet]=Update(state,inlet,outlet,dt,zeta)
%---------------------Update Extrinsic Properties-------------------------%
state.enthalpy=state.enthalpy+state.heataddition;
state.heataddition(:)=0;

switch state.type
    case '2phasefluid'
        %--------------------------Mass-Exchange--------------------------%
        [state,outlet]=transport(outlet,inlet,state,dt,zeta);
        %----------------------Calculate Pressure-------------------------%
        state=Pressure9(state);
        %--------------------Update Intrinsic Properties------------------%
        state=Props(state);
        %--------------------Update Element Distribution------------------%
        state=remesh(state);
        %---------------------Calculate Values at Outet-------------------%
        outlet.fluid=state.fluid;
        outlet.type=state.type;
        outlet.dimensions=0;
        outlet=Props(outlet);
    case '1phasefluid'
        %--------------------------Mass-Exchange--------------------------%
        [state,outlet]=transport(outlet,inlet,state,dt,zeta);
        %----------------Update Intrinsic Properties----------------------%
        state=Props(state);
        %---------------------Calculate Values at Outet-------------------%
        outlet.fluid=state.fluid;
        outlet.type=state.type;
        outlet.dimensions=0;
        outlet=Props(outlet);
    case 'solid'
        M_w=(state.mass.*state.Cp);
        %----------------Update Intrinsic Properties----------------------%
        state.temperature=state.enthalpy./M_w;
        %---------------------Calculate Values at Outet-------------------%
        outlet=[];
end
