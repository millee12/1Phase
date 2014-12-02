%Single-phase Case
clear all
tic
clc
warning('off','MATLAB:interp1:NaNinY')
set(0,'defaulttextinterpreter','latex')
%#ok<*SAGROW>
%% -----------------------Simulation Parameters---------------------------%
%time
global zeta dt t
dt = 0.05;
time=0:dt:200;
%remeshing tolerance
zeta=.001;
%% --------------Component Geometries, Initial Conditions-----------------%
%Initial conditions
T_amb=75+460;
%=====Cold Fluid=====%
state_C.name='Cold Fluid';
state_C.type='1phasefluid';
state_C.dimensions=1;
state_C.N=50;
%Geometry
state_C.entrance=1;
state_C.exit=state_C.N;
state_C.geometry.breakpoints=[0,50];
state_C.geometry.crosssections=.5;
state_C=firstmesh(state_C,0,'PAO',T_amb);
%Properties
state_C.resistance=5E3*ones(state_C.N,1);
%=====Hot fluid=======%
state_H.name='Hot Fluid';
state_H.type='1phasefluid';
state_H.dimensions=1;
state_H.N=50;
%Geometry
state_H.entrance=state_H.N;
state_H.exit=1;
state_H.geometry.breakpoints=[0,50];
state_H.geometry.crosssections=2.5;
state_H.resistance=2E3*ones(state_H.N,1);
state_H=firstmesh(state_H,0,'PAO',T_amb);
%=====Heat Exchanger=======%
state_HX.N=25;
tmp.name='Aluminum2024';
R=load('fluids.mat',tmp.name);
tmp=R.(genvarname(tmp.name));
state_HX.type='solid';
state_HX.dimensions=1;
state_HX.mass=(3/state_HX.N)*ones(state_HX.N,1);
state_HX.Cp=tmp.Cp*ones(state_HX.N,1);
state_HX.enthalpy=T_amb.*state_HX.mass.*state_HX.Cp;
state_HX.temperature=T_amb*ones(state_HX.N,1);
state_HX.resistance=.25/tmp.k*ones(state_HX.N,1);
state_HX.inside.perimeter=20;
state_HX.outside.perimeter=20;
state_HX.mesh=linspace(0,50,state_HX.N+1)';
TmeanHX=zeros(size(time));
%% ------------------------Boundary Conditions----------------------------%
%Hot Fluid
inlet_H.temperature=100+460;
fluid=state_H.fluid;
inlet_H.massflow=.1;
inlet_H.specificenthalpy=interp1DUG(fluid.TT,fluid.H_T,inlet_H.temperature);
outlet_H.massflow=inlet_H.massflow;
ToutH=zeros(size(time));
%Cold Fluid
inlet_C.temperature=50+460;
fluid=state_C.fluid;
inlet_C.massflow=.05;
inlet_C.specificenthalpy=interp1DUG(fluid.TT,fluid.H_T,inlet_C.temperature);
outlet_C.massflow=inlet_C.massflow;
ToutC=zeros(size(time));
%% ----------------------------Simulate-----------------------------------%
i=1;
%---------------------------------Model-----------------------------------%
tic
for t=time
    if mod(t,time(end)/100)==0
        %Counter
        dispstat(['Simulation ',...
            num2str(100*(t-time(1))/(time(end)-time(1))),'% Complete'])
    end
    % ______________________Heat Exchange________________________________%
    %++++++++++++++++++++++++++LOW SIDE+++++++++++++++++++++++++++++++++++%
    %Hot Fluid <---> HX
    [state_H,state_HX,QdotH]=HX2(state_H,state_HX,'inside');
    %Cold Fluid <---> HX
    [state_C,state_HX,QdotC]=HX2(state_C,state_HX,'outside');
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
    % ______________________Mass Exchange________________________________%
    [state_H,outlet_H]=transport(outlet_H,inlet_H,state_H);
    [state_C,outlet_C]=transport(outlet_C,inlet_C,state_C);
    % --------------------Update Intrinsic Properties------------------%
    state_H=Props(state_H);
    state_C=Props(state_C);
    state_HX=Props(state_HX);
    % --------------------Update Element Distribution------------------%
    state_H=remesh(state_H);
    state_C=remesh(state_C);
    % _____________________Store Information_____________________________%
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
% plot(xHX,THX,'k',xC,TC,'b','LineWidth',2)
xlim([xHX(1),xHX(end)])
% legend('fluid A','Heat Exchanger','fluid B','Location','Southwest')
legend('Hot Fluid','Heat Exchanger','Cold Fluid','Location','Southwest')
xlabel('x (in)')
ylabel('temperature (\(^{o}\) F)')
legend('boxoff')
saveas(gcf,'fluidsx','epsc')

figure(2);
set(gca,'Fontsize',20,'LineWidth',2)
plot(time,ToutH,'r',time,TmeanHX,'k',time,ToutC,'b','LineWidth',2)
% plot(time,TmeanHX,'k',time,ToutC,'b','LineWidth',2)
xlim([0,time(end)])
xlabel('time (s)')
ylabel('temperature (\(^{o}\) F)')
legend('fluid A outlet','Heat Exchanger mean','fluid B outlet','Location','Southeast')
% legend('Solid, mean','Fluid, outlet','Location','northeast')
legend('boxoff')
saveas(gcf,'fluidst','epsc')

figure(3);
C=jet;
C(1,:)=1;
set(gca,'Fontsize',20,'LineWidth',2)
imagesc(QdotH)
colormap(C)
h=colorbar;
set(h,'fontsize',14);
xlabel('Heat Exchanger')
ylabel('Hot Fluid')
saveas(gcf,'qdoth','epsc')

figure(4);
C=jet;
C(end,:)=1;
set(gca,'Fontsize',20,'LineWidth',2)
imagesc(QdotC)
xlim([1,state_HX.N])
ylim([1,state_C.N])
colormap(C)
h=colorbar;
set(h,'fontsize',14);
xlabel('Solid (b)')
ylabel('Cold Fluid (a)')
title(['t = ',  num2str(t),'s','\(\dot{Q}_{a \to b}\) (BTU/s) '])
saveas(gcf,'qdotc','epsc')

tt=toc;