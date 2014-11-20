%% Initial Conditions For Full VCS cycle
%% Model Agents
%The model may involve any number of solid or liquid materials.  These will
%be defined as agents and must be given spatial bounds on the 1-D line that
%represents the entire pressure side (high or low)
%Finite Volume Lengths
% Set to zero to remove pressure side from simulation

%mesh size
N_H=20;
N_C=20;
N_HX=100;
%ambient conditions
T_amb=75+460;
%Solids
start_HX=0;
L_HX=50;
%Fluids
start_H=0;
L_H=50;
start_C=0;
L_C=50;
%% Exchanger Wall
wall_mat='Aluminum2024';
R=load('fluids.mat',wall_mat);
wall_props=R.(genvarname(wall_mat));
wall_props.name=wall_mat;
clear wall_mat

k_HX=wall_props.k;       %BTU/(s in R)
rho_HX=wall_props.rho;   %lb/in3
Cp_HX=wall_props.Cp;     %BTU/(lb R)

%Passes
%wall mass/volume
m_HX=1/N_HX;                          %lbm/in (mass/in)
%inner Surface Area
per_HX_i=20;                       %in^2/l (surface area per length)
%Outer Surface Area
per_HX_o=20;
%Cross-Sectional flow Area
Acs_HX_i=5;                         %in^2
Acs_HX_o=5;                           %in^2
%state_HX=zeros(N_HX,9);
T_0_HX=T_amb*ones(length(N_HX),1);
state_HX.type='solid';
state_HX.mass=m_HX*ones(N_HX,1);
state_HX.length=(L_HX/N_HX)*ones(N_HX,1);
state_HX.enthalpy=T_0_HX.*m_HX.*Cp_HX*ones(N_HX,1);
state_HX.temperature=T_0_HX*ones(N_HX,1);
state_HX.density=rho_HX*ones(N_HX,1);
state_HX.heataddition=zeros(N_HX,1);
state_HX.alpha=k_HX*ones(N_HX,1);
state_HX.Cp=Cp_HX*ones(N_HX,1);
state_HX.inside.start=start_HX;
state_HX.inside.totallength=L_HX;
state_HX.inside.perimeter=per_HX_i;
state_HX.outside.start=start_HX;
state_HX.outside.totallength=L_HX;
state_HX.outside.perimeter=per_HX_o;
state_HX.inside.thickness=.25;
state_HX.outside.thickness=.25;
state_HX.mesh=linspace(start_HX,start_HX+L_HX,N_HX+1)';
state_HX.N=N_HX;
disp(['Heat Exchanger Loaded...', wall_props.name,])


%% Generate Fluid Meshes
%Hot fluid
state_H.geometry.breakpoints=[start_H,L_H];
state_H.geometry.crosssections=Acs_HX_i;
state_H.N=N_H;
state_H.type='1phasefluid';
state_H.entrance=1;
state_H.exit=N_H;
state_H.dimensions=1;
state_H.name='Hot Fluid';
state_H.alpha=.001*ones(N_H,1);
state_H=firstmesh(state_H,0,'PAO',T_amb);


disp(['Hot Fluid Loaded...', state_H.fluid.name])
%Cold Fluid
state_C.geometry.breakpoints=[start_C,L_C];
state_C.geometry.crosssections=Acs_HX_o;
state_C.N=N_C;
state_C.type='1phasefluid';
state_C.entrance=N_C;
state_C.exit=1;
state_C.dimensions=1;
state_C.name='Cold Fluid';
state_C.alpha=.001*ones(N_C,1);
state_C=firstmesh(state_C,0,'PAO',T_amb);

disp(['Cold Fluid Loaded...', state_C.fluid.name])



