%Single-phase Case
clear all
tic
clc
warning('off','MATLAB:interp1:NaNinY')
set(0,'defaulttextinterpreter','latex')
%#ok<*SAGROW>
%% -----------------------Simulation Parameters---------------------------%
%time
global dt t
dt = 0.05;
time=0:dt:125;
%remeshing tolerance
zeta=.001;
%% --------------Component Geometries, Initial Conditions-----------------%
%Initial conditions
T_amb=75+460;
%% Cold Fluid
state_C.N=75;
%Fluid Type
state_C.type='1phasefluid';
%Flow Direction
state_C.entrance=state_C.N;
state_C.exit=1;
%Cold Fluid
state_C.geometry.breakpoints=[0,50];
state_C.geometry.crosssections=5;
state_C.dimensions=1;
state_C.name='Cold Fluid';
state_C.resistance=5000*ones(state_C.N,1);
state_C=firstmesh(state_C,0,'PAO',T_amb);
%% Hot fluid
state_H.N=150;
state_H.type='1phasefluid';
%Flow Direction
state_H.entrance=1;
state_H.exit=state_H.N;
%Geometry
state_H.geometry.breakpoints=[0,50];
state_H.geometry.crosssections=5;
state_H.dimensions=1;
state_H.name='Hot Fluid';
state_H.resistance=2000*ones(state_H.N,1);
state_H=firstmesh(state_H,0,'PAO',T_amb);
%% Heat Exchanger
state_HX.N=100;
tmp.name='Aluminum2024';
R=load('fluids.mat',tmp.name);
tmp=R.(genvarname(tmp.name));
state_HX.type='solid';
state_HX.dimensions=1;
state_HX.mass=(1.5/state_HX.N)*ones(state_HX.N,1);
state_HX.Cp=tmp.Cp*ones(state_HX.N,1);
state_HX.enthalpy=T_amb.*state_HX.mass.*state_HX.Cp;
state_HX.temperature=T_amb*ones(state_HX.N,1);
state_HX.resistance=.25/tmp.k*ones(state_HX.N,1);
state_HX.inside.perimeter=20;
state_HX.outside.perimeter=20;
state_HX.mesh=linspace(5,45,state_HX.N+1)';
%% ------------------------Boundary Conditions----------------------------%
%Hot Fluid
inlet_H.temperature=100+460;
fluid=state_H.fluid;
inlet_H.massflow=.5;
inlet_H.specificenthalpy=interp1DUG(fluid.TT,fluid.H_T,inlet_H.temperature);
outlet_H.massflow=inlet_H.massflow;
%Cold Fluid
inlet_C.temperature=50+460;
fluid=state_C.fluid;
inlet_C.massflow=.1;
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
    if mod(t,time(end)/100)==0
        %Counter
        dispstat(['Simulation ',...
            num2str(100*(t-time(1))/(time(end)-time(1))),'% Complete'])
    end
    %% ______________________Heat Exchange________________________________%
    %++++++++++++++++++++++++++LOW SIDE+++++++++++++++++++++++++++++++++++%
    %Hot Fluid <---> HX
    [state_H,state_HX,QdotH]=HX2(state_H,state_HX,'inside');
    %Cold Fluid <---> HX
    [state_C,state_HX,QdotC]=HX2(state_C,state_HX,'outside');
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
    %% ______________________Mass Exchange________________________________%
    [state_H,outlet_H]=transport(outlet_H,inlet_H,state_H,dt,zeta);
    [state_C,outlet_C]=transport(outlet_C,inlet_C,state_C,dt,zeta);
    %% --------------------Update Intrinsic Properties------------------%
    state_H=Props(state_H);
    state_C=Props(state_C);
    state_HX=Props(state_HX);
    %% --------------------Update Element Distribution------------------%
    state_H=remesh(state_H);
    state_C=remesh(state_C);
    %% _____________________Store Information_____________________________% 
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
xlim([xHX(1),xHX(end)])
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
legend('fluid A outlet','Heat Exchanger mean','fluid B outlet','Location','Southeast')
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
colormap(C)
h=colorbar;
set(h,'fontsize',14);
xlabel('Heat Exchanger')
ylabel('Cold Fluid')
saveas(gcf,'qdotc','epsc')

tt=toc;