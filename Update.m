function [state]=Update(state,dt)
%---------------------Update Extrinsic Properties-------------------------%
state.enthalpy=state.enthalpy+state.heatrate*dt;
state.heatrate(:)=0;
%----------------------Calculate Pressure-------------------------%
state=Pressure9(state);
%--------------------Update Intrinsic Properties------------------%
state=Props(state);
%--------------------Update Element Distribution------------------%
state=remesh(state);

