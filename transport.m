function [state,outlet]=transport(outlet,inlet,state,dt,zeta)
%% Assign Variables
n=state.N;
%inlet
mdot_in=inlet.massflow;
h_in=inlet.specificenthalpy;
%outlet
mdot_out=outlet.massflow;
%fluid properties
h=state.enthalpy./state.mass;
entrance=state.entrance;
exit=state.exit;
nee=length(exit);
%Count amalgamations
%% Empty Exit Volume
outlet.massflow=mdot_out;
state.mass(exit)=state.mass(exit)-(mdot_out/nee)*dt;
outlet.specificenthalpy=mean(h(exit));
state.enthalpy(exit)=...
    state.enthalpy(exit)-outlet.specificenthalpy*(mdot_out/nee)*dt;
%% Fill Entrance Volume
state.enthalpy(entrance)=state.enthalpy(entrance)+(mdot_in/nee)*h_in*dt;
state.mass(entrance)=state.mass(entrance)+(mdot_in/nee)*dt;
%% Remove then add Elements
minmass=(sum(state.mass)/n)*zeta;
w=1;
wmax=n;
mn=sum(state.mass);
while min(state.mass) < minmass
    %Find Smallest Element 
    i=find(state.mass < minmass,1);
    %Add its mass and enthalpy to adjoining elements
    merge=[i-1;i+1];
    if merge(1) < 1
        merge(1)=i+2;
    elseif merge(2) > n
        merge(2)=i-2;
    end
    state.mass(merge)=state.mass(merge)+state.mass(i)/2;
    state.enthalpy(merge)=state.enthalpy(merge)+state.enthalpy(i)/2;
    %Remove Smallest Element
    state.mass=[state.mass(1:i-1);state.mass(i+1:n)];
    state.enthalpy=[state.enthalpy(1:i-1);state.enthalpy(i+1:n)];
    %Find Largest Element
    j=find(state.mass==max(state.mass),1);
    %Split it to create a new element
    state.mass(j)=state.mass(j)/2;
    state.enthalpy(j)=state.enthalpy(j)/2;
    
    state.mass=[state.mass(1:j);state.mass(j:n-1)];
    state.enthalpy=[state.enthalpy(1:j);state.enthalpy(j:n-1)]; 
    if w > wmax
        error(['failure in transport for ', state.name])
    end
    w=w+1;
end
%Error Checking
if length(state.mass) ~= n
    error(['transport changed element number in',state.name])
elseif sum(state.mass) <= mn -0.01
    error(['transport changed mass in ',state.name])
end
%---------------------Calculate Values at Outlet-------------------%
outlet.fluid=state.fluid;
outlet.type=state.type;
outlet.dimensions=0;
outlet=Props(outlet);

