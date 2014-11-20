%Single-phase Case
clear all
tic
clc
warning('off','MATLAB:interp1:NaNinY')
set(0,'defaulttextinterpreter','latex')
%#ok<*SAGROW>
%% -----------------------Simulation Parameters---------------------------%
%time
dt=0.01;
time=0:dt:50;
global t
%remeshing tolerance
zeta=.001;
%% --------------Component Geometries, Initial Conditions-----------------%
geometry
%% ------------------------Boundary Conditions----------------------------%
%Hot Fluid
inlet_H.temperature=90+460;
fluid=state_H.fluid;
inlet_H.massflow=1;
inlet_H.specificenthalpy=interp1DUG(fluid.TT,fluid.H_T,inlet_H.temperature);
outlet_H.massflow=inlet_H.massflow;
%Cold Fluid
inlet_C.temperature=50+460;
fluid=state_C.fluid;
inlet_C.massflow=.25;
inlet_C.specificenthalpy=interp1DUG(fluid.TT,fluid.H_T,inlet_C.temperature);
outlet_C.massflow=inlet_C.massflow;
%% ----------------------------Simulate-----------------------------------%
i=1;
%---------------------------------Model-----------------------------------%
tic
for t=time
    if mod(t,time(end)/10)==0
        %Counter
        disp(['Simulation ',...
            num2str(100*(t-time(1))/(time(end)-time(1))),'% Complete'])
    end
    %% ______________________Heat Exchange________________________________%
    %++++++++++++++++++++++++++LOW SIDE+++++++++++++++++++++++++++++++++++%
    %Hot Fluid <---> HX
    [state_H,state_HX]=HX(state_H,state_HX,'inside',dt);
    %Cold Fluid <---> HX
    [state_C,state_HX]=HX(state_C,state_HX,'outside',dt);
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%___________________________%
    %% ______________________Update Sections_____________________________%
    %Hot fluid
    [state_H,outlet_H]=Update(state_H,inlet_H,outlet_H,dt,zeta);
    %solid
    [state_HX,~]=Update(state_HX,0,0,dt,0);
    %Cold fluid
    [state_C,outlet_C]=Update(state_C,inlet_C,outlet_C,dt,zeta);
    i=i+1;
end
%% Plot
xH=[state_H.geometry.breakpoints(1);...
    (state_H.mesh(1:end-1)+state_H.mesh(2:end))/2;...
    state_H.geometry.breakpoints(end)];
xC=[state_C.geometry.breakpoints(1);...
    (state_C.mesh(1:end-1)+state_C.mesh(2:end))/2;...
    state_C.geometry.breakpoints(end)];

TH=[inlet_H.temperature;state_H.temperature;outlet_H.temperature]-460;
TC=[outlet_C.temperature;state_C.temperature;inlet_C.temperature]-460;

figure(1);plot(xC,TC,xH,TH)
i=i-1;
tt=toc;