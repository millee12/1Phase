%% Interior fluid-Low Side
if L_i_LS>0
    %Total Length
    L=L_i_LS;
    %Allocate initial matrix
    initial_i_LS=zeros(N_i_LS,9);
    % Total Volume
    V_i_LS=(L_p1)*Acs_p1+(L_evap)*Acs_evap_i+(L_p2)*Acs_p2;
    v_i_LS=V_i_LS/(L_p1+L_evap+L_p2);
    %Choose Interior Fluid
    INT_FLUID='R134a';
    R=load('fluids.mat',INT_FLUID);
    fluid_i=R.(genvarname(INT_FLUID));
    fluid_i.name=INT_FLUID;
    %1.  Masses
    m_i_LS=m_0_LS/(N_i_LS);
    M_n=m_i_LS*ones(N_i_LS,1);
    %2.  Volume Lengths
    l_n=L/(N_i_LS-1);
    l_n_ends=l_n/2;
    L_n=l_n*ones(N_i_LS,1);
    %3.  Pressure
    %(saturation pressure for T_amb)
    P=interp2(fluid_i.TT,fluid_i.QQ,fluid_i.P_T_Q,T_amb,.5);
    T_sat=interp2(fluid_i.PP,fluid_i.QQ,fluid_i.T_P_Q,P,.5);
    %Saturated Values
    h_v=interp2(fluid_i.PP,fluid_i.QQ,fluid_i.U_P_Q,P,.99);
    h_l=interp2(fluid_i.PP,fluid_i.QQ,fluid_i.U_P_Q,P,0);
    rho_v=interp2(fluid_i.PP,fluid_i.QQ,fluid_i.D_P_Q,P,.99);
    rho_l=interp2(fluid_i.PP,fluid_i.QQ,fluid_i.D_P_Q,P,0);
    %7.  Density
    rho=(m_0_LS/V_i_LS);
    %6.  Quality
    Q=((rho^-1)-(rho_l^-1))/((rho_v^-1)-(rho_l^-1));
    %4.  Enthalpies
    %specific
    h_i_LS=h_l+Q*(h_v-h_l);
    %Total
    H=h_i_LS.*m_i_LS*ones(N_i_LS,1);
    %Fluid Matrix
    initial_i_LS(:,:)=[M_n L_n ones(N_i_LS,1)*P H ones(N_i_LS,1)*[T_amb Q rho 0 0]];
    initial_i_LS=props2P(fluid_i,initial_i_LS,'states');
    %Geometry Vector
    props_i_LS=[L_i_LS 0 v_i_LS start_i_LS 0 0];
    geo_i_LS=[pipe1_props;evap_props(1,:);pipe2_props];
    %Feedback
    disp(['Low-Side Interior Fluid Loaded...', fluid_i.name,])
end
%% Exterior fluid-Low Side
if L_o_LS>0
    EXT_FLUID='Air';
    R=load('fluids.mat',EXT_FLUID);
    fluid_o_LS=R.(genvarname(EXT_FLUID));
    fluid_o_LS.name=EXT_FLUID;
    clear EXT_FLUID
    
    initial_o_LS=zeros(N_o_LS,9);
    P_0=11;
    % Total Volume
    V_o_LS=L_evap*Acs_evap_o;
    v_o_LS=V_o_LS/N_o_LS;
    %Volume Lengths
    initial_o_LS(:,2)=(L_o_LS/N_o_LS)*ones(N_o_LS,1);
    %Pressures
    initial_o_LS(:,3)=P_0*ones(N_o_LS,1);
    %Masses
    initial_o_LS(:,1)=initial_o_LS(:,2).*v_o_LS.*...
        interp1(fluid_o_LS.TT,fluid_o_LS.D_T,T_amb);
    %Enthalpies
    initial_o_LS(:,4)=initial_o_LS(:,1).*...
        interp1(fluid_o_LS.TT,fluid_o_LS.H_T,T_amb);
    %Intrinsic Values
    [initial_o_LS]=props1P(fluid_o_LS,initial_o_LS,'states');
    props_o_LS=[L_o_LS 0 v_o_LS start_o_LS];
    geo_o_LS=evap_props(2,:);
    geo_o_LS(1)=L_o_LS;
    disp(['Low-Side Exterior Fluid Loaded...', fluid_o_LS.name,])
else
    initial_o_LS=[];
    state_o=[];
end
%% Interior fluid-High Side
if L_i_HS>0
    INT_FLUID='R134a';
    R=load('fluids.mat',INT_FLUID);
    fluid_i=R.(genvarname(INT_FLUID));
    fluid_i.name=INT_FLUID;
    V_i_HS=(L_p3)*Acs_p3+(L_cond)*Acs_cond_i+(L_p4)*Acs_p4;
    v_i_HS=V_i_HS/(L_p3+L_cond+L_p4);
    %1.  Mass
    m_i_HS=(m_0_HS/N_i_HS)*ones(N_i_HS,1);
    %2. Length
    l_n=L_i_HS/N_i_HS*ones(N_i_HS,1);
    %3.  Pressure
    %(saturation pressure for T_amb)
    P=interp2DUG(fluid_i.TT,fluid_i.QQ,fluid_i.P_T_Q,T_amb,.5);
    %Saturated Values
    sat=interp1DUG(fluid_i.PP,fluid_i.sat_par,P');
    h_l=sat(:,1);
    h_v=sat(:,2);
    rho_l=sat(:,3);
    rho_v=sat(:,4);
    %7.  Density
    rho=(m_0_HS./V_i_HS)*ones(N_i_HS,1);
    %6.  Quality
    Q=((rho.^-1)-(rho_l.^-1))./((rho_v.^-1)-(rho_l.^-1));
    %4.  Enthalpies
    %specific
    h_i_HS=h_l+Q*(h_v-h_l);
    %Total
    H=h_i_HS.*m_i_HS;
    %Fluid Matrix
    initial_i_HS=[m_i_HS l_n ones(N_i_HS,1)*P H ones(N_i_HS,1)*T_amb Q rho zeros(N_i_HS,1) zeros(N_i_HS,1)];
    %Geometry Vector
    props_i_HS=[L_i_HS 0 v_i_HS start_i_HS];
    geo_i_HS=[pipe3_props;cond_props(1,:);pipe4_props];
    %Feedback
    disp(['High-Side Interior Fluid Loaded...', fluid_i.name,])
end
%% Exterior fluid-High Side
if L_o_HS>0
    V_o_HS=L_cond*Acs_cond_o;
    v_o_HS=V_o_HS/N_o_HS;
    l_n=(L_o_HS/N_o_HS)*ones(N_o_HS,1);
    %Choose Interior Fluid
    INT_FLUID='Air';
    R=load('fluids.mat',INT_FLUID);
    fluid_o_HS=R.(genvarname(INT_FLUID));
    fluid_o_HS.name=INT_FLUID;
    %3.  Pressure
    P_0=0;
    %7.  Density
    rho=interp1(fluid_o_HS.TT,fluid_o_HS.D_T,T_amb)*ones(N_o_HS,1);
    %1.  Mass
    m_o_HS=(rho*v_o_HS);
    %6.  Quality
    Q=1*ones(N_o_HS,1);
    %4.  Enthalpies
    %specific
    h_o_LS=interp1(fluid_o_HS.TT,fluid_o_HS.H_T,T_amb)*ones(N_o_HS,1);
    %Total
    H=h_o_LS.*m_o_HS;
    %Geomerty Vector
    props_o_HS=[L_o_HS 0 v_o_HS start_o_HS ];
    geo_o_HS=cond_props(2,:);
    geo_o_HS(1)=L_o_HS;
    %Fluid Matrix
    initial_o_HS=[m_o_HS l_n ones(N_o_HS,1)*P_0 H ones(N_o_HS,1)*T_amb Q rho zeros(N_o_HS,1) zeros(N_o_HS,1)];
    %Feedback
    disp(['High-Side Exterior Fluid Loaded...', fluid_o_HS.name,])
end
%% Load Steady State Condtions
if steady==1
    [initial_i_LS,initial_i_HS, initial_o_LS,initial_o_HS,...
        initial_w_LS, initial_w_HS, pp_init, rpm_init, mdot_water_init]=...
        steadystart('./stead1');
else
    pp_init=0;
    rpm_init=0;
    mdot_water_init=0;
end
z=length(time)+1;