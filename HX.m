%Heat Exchanger Model
function [state_a,state_b] = HX(state_a,state_b,location,dt)
%% Define Variables
%state a
if strcmp(state_a.type,'solid')
    per=eval(['state_a.',location,'.perimeter']);
    lambda_a=eval(['state_a.',location,'.thickness']);
else
    lambda_a=1;
end
%state b
if strcmp(state_b.type,'solid')
    per=eval(['state_b.',location,'.perimeter']);
    lambda_b=eval(['state_b.',location,'.thickness']);
else
    lambda_b=1;
end
%% Calculate Interfaces
[l_int]=calcinter2(state_a.mesh,state_b.mesh);
%% Calculate Thermal Resistance
U_a=ones(state_b.N,1)*state_a.alpha';
U_b=state_b.alpha*ones(1,state_a.N);
As=l_int.*per;
%% Generate Temperature Gradient Matrix
T_a=state_a.temperature;                
T_b=state_b.temperature;
DT=T_b*ones(1,state_a.N)-ones(state_b.N,1)*T_a';
R=(lambda_b./(As.*U_b))+(lambda_a./(As.*U_a));
Qdot=DT./R;
%% Calculate Heat-Transfer Rate to each finite-volume
dH_a=sum(Qdot,1)*dt;         % Enthalpy flowrate for Each Volume
dH_b=-sum(Qdot,2)'*dt;       % Energy flowrate for Each Volume
%% Assemble Output Matrices
state_a.heataddition=state_a.heataddition+dH_a';
state_b.heataddition=state_b.heataddition+dH_b';
    