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
time=0:dt:150;
global t
%remeshing tolerance
zeta=.001;
%% --------------Component Geometries, Initial Conditions-----------------%
%mesh size
N_H=50;
N_C=50;
N_HX=25;
%Initial conditions
T_amb=75+460;
T0_HX=50+460;
T0_H=100+460;
T0_C=50+460;
%Fluid Type
state_C.type='1phasefluid';
%Flow Direction
state_C.entrance=N_C;
state_C.exit=1;
%Flow Direction
state_H.type='1phasefluid';
state_H.entrance=1;
state_H.exit=N_H;
%Run Geometry Script
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
%% Initialize
    ToutH=zeros(size(time));
    TmeanHX=zeros(size(time));
    ToutC=zeros(size(time));
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
    ToutH(i)=outlet_H.temperature-460;
    TmeanHX(i)=mean(state_HX.temperature)-460;
    ToutC(i)=outlet_C.temperature-460;
    i=i+1;
end
%% Plotting
xH=(state_H.mesh(1:end-1)+state_H.mesh(2:end))/2;
xHX=(state_HX.mesh(1:end-1)+state_HX.mesh(2:end))/2;
xC=(state_C.mesh(1:end-1)+state_C.mesh(2:end))/2;

TH=state_H.temperature-460;
THX=state_HX.temperature-460;
TC=state_C.temperature-460;

figure(1);
set(gca,'Fontsize',20,'LineWidth',2)
plot(xH,TH,'r',xHX,THX,'k',xC,TC,'b','LineWidth',2)
legend('fluid A','Heat Exchanger','fluid B','Location','Southwest')
xlabel('x (in)')
ylabel('temperature (\(^{o}\) F)')
legend('boxoff')
saveas(gcf,'fluidsx','epsc')

figure(2);
set(gca,'Fontsize',20,'LineWidth',2)
plot(time,ToutH,'r',time,TmeanHX,'k',time,ToutC,'b','LineWidth',2)
xlabel('time (s)')
ylabel('temperature (\(^{o}\) F)')
legend('fluid A outlet','Heat Exchanger mean','fluid B outlet','Location','Southwest')
legend('boxoff')
saveas(gcf,'fluidst','epsc')
tt=toc;