function [acc_z] = acc_func(X)
    %% Add Path
    % Current dir
    currentDir = fileparts(mfilename('fullpath'));

    % Add the 'Util' path
    utilPath = fullfile(currentDir, 'Util');
    addpath(utilPath);
    %% Handles Input
    percent_foam = X(1);
    percent_metal = X(2);

    %% Initial Speed and Position in NED Frame
    IC_Pos = [0; 0; 10; 0; 0; 0];
    IC_Velo = [0; 0; 0; 0; 0; 0];

    %% General Parameters
    % Environment
    rho_fluid = 1000;      % Fresh water density, kg/m3
    mu_fluid = 1.002e-3;   % Fresh water dynamic viscosity, Pa.s
    g = 9.81;              % m/s2, Positive pointing down to center of Earth

    % Sphere Properties
    D = 2;                                     % Diameter, m
    R = D/2;                                % Radius, m
    V = 4/3 * pi * R^3;                     % m3

    rho_foam = 300;                               % Heavy Foam, kg/m3
    rho_metal = 7500;                             % Stainless Steel, kg/m3
    rho_cavity = rho_fluid;                             % kg/m3

    foam_height = percent_foam * R;
    metal_height = percent_metal * R;
    c_t_height = R - foam_height;
    c_b_height = R - metal_height;

    V_foam = pi * (foam_height^2) *(R - foam_height/3);
    V_metal = pi * (metal_height^2) *(R - metal_height/3);

    % V_foam = pi*((R^3 - R * c_t_height^2) - ((R^3 - c_t_height^3)/3));
    % V_metal = pi*((R^3 - R * c_b_height^2) - ((R^3 - c_b_height^3)/3));

    V_c_t = (V/2) - V_foam;
    V_c_b = (V/2) - V_metal;
    V_cavity = V_c_t + V_c_b;

    m_foam = rho_foam * V_foam;        % kg
    m_metal = rho_metal * V_metal;     % kg
    m_cavity = rho_cavity * V_cavity;  % kg

    m = m_foam + m_metal + m_cavity;   % kg
    mb = rho_fluid * V;                       % kg
    I = 2/5 * m * R^2;                        % Moment inertia, kg.m2
    ma = -rho_fluid * 2/3 * pi * R^3;         % Added mass, kg

    W = m * g;                    % N
    B = mb * g;                   % Fully submerged, N

    %% Center of Origin
    % Center of Origin defined in Body Frame
    xo = 0;
    yo = 0;
    zo = 0;
    
    ro = [xo; yo; zo];

    %% Centroid in z-axis for each Parts
    % Foam
    chord_f = 2*sqrt(foam_height*(2*R-foam_height));
    gamma_f = 2*asin(chord_f/(2*R));
    num_foam = 4*R*(sin(gamma_f/2))^3;
    denum_foam = 3*(gamma_f - sin(gamma_f));
    if foam_height ~= 0
        centroid_foam = num_foam / denum_foam;
    else
        centroid_foam = 0;
    end

    %centroid_foam = ((((R^4)/2 - (R^2 * c_t_height^2)/2) - ((R^4)/4 - (c_t_height^4)/4))/V_foam + c_t_height);

    % Metal
    chord_m = 2*sqrt(metal_height*(2*R-metal_height));
    gamma_m = 2*asin(chord_m/(2*R));
    num_metal = 4*R*(sin(gamma_m/2))^3;
    denum_metal = 3*(gamma_m - sin(gamma_m));
    if metal_height ~= 0
        centroid_metal = -num_metal / denum_metal;
    else
        centroid_metal = 0;
    end

    %centroid_metal = -((((R^4)/2 - (R^2 * c_b_height^2)/2) - ((R^4)/4 - (c_b_height^4)/4))/V_foam + c_b_height);

    % Cavity
    centroid_HS = 3 * R / 8;
    V_HS = V/2;

    centroid_cavity_t = ((centroid_HS * V_HS) - (centroid_foam * V_foam))/V_c_t;
    centroid_cavity_b = -((centroid_HS * V_HS) - (centroid_metal * V_metal))/V_c_b;
    centroid_cavity = ((centroid_cavity_t * V_c_t) + (centroid_cavity_b * V_c_b))/(V_cavity);

    %% Center of Gravity and Buoyancy in Body Frame
    % Center of Gravity
    xg = 0;
    yg = 0;
    zg = (m_foam*centroid_foam + m_metal*centroid_metal + m_cavity*centroid_cavity)/(m_foam + m_metal + m_cavity);
    
    rg = [xg; yg; zg];

    % Center of Buoyancy
    xb = 0;
    yb = 0;
    zb = (V_foam*centroid_foam + V_metal*centroid_metal + V_cavity*centroid_cavity)/(V);

    rb = [xb; yb; zb];

    % Eulerian distance for transformation
    rg_o = rg - ro; % Distance from center of gravity to center of origin 
    rb_o = rb - ro; % Distance from center of buoyancy to center of origin

    %% Body Positions
    % Only used if body is modular and have several components

    %% Body Mass Matrix Transformed to CO
    Mrb = [m 0 0 0 m*zg 0;
             0 m 0 -m*zg 0 0;
             0 0 m 0 0 0;
             0 -m*zg 0 I 0 0;
             m*zg 0 0 0 I 0;
             0 0 0 0 0 I];

    Mrb_o = Transform(Mrb,rg_o);

    %% Body Added Mass Matrix Transformed to CO
    Ma = -diag([ma ma ma 0 0 0]);

    Ma_o = Transform(Ma,rb_o);

    %% Generalized Mass Matrix
    MT = Mrb_o + Ma_o;
    
    %% Position
    % Described in NED frame
    x = IC_Pos(1,1);
    y = IC_Pos(2,1);
    z = IC_Pos(3,1);
    phi = IC_Pos(4,1);
    theta = IC_Pos(5,1);
    psi = IC_Pos(6,1);

    % General Vector
    T_ = [x; y; z];         % Translational Vector
    R_ = [phi, theta, psi]; % Rotational Vector

    %% Velocity
    % Described in Origin frame
    u = IC_Velo(1,1);
    v = IC_Velo(2,1);
    w = IC_Velo(3,1);
    p = IC_Velo(4,1);
    q = IC_Velo(5,1);
    r = IC_Velo(6,1);

    % General Vector
    V_ = [u; v; w];     % Translation Velocity vector, read "Nu"
    W_ = [p; q; r];     % Angular velocity vector, read "Omega"

    V_t = [V_ ; W_];    % 6DOF Velocity

    % %% Gravity Force
    % %Described in Origin Frame
    % Fg_F = [W * sin(theta);
    %         -W * cos(theta) * sin(phi);
    %         -W * cos(theta) * cos(phi)];  % Gravity Forces
    % Fg_M = zeros(3,1);                          % Gravity Moment
    % 
    % Fg = [Fg_F ; Fg_M];
    % 
    % %% Buoyancy Force
    % %Described in Body Frame (Depends on the submerged volume)
    % Fb_F = [-B * sin(theta);
    %         B * cos(theta) * sin(phi);
    %         B * cos(theta) * cos(phi)];   % Buoyancy Forces
    % lever_arm = -rb + rg;            
    % Fb_M = SS_(lever_arm) * Fb_F;               % Buoyancy Moment
    % 
    % Fb = [Fb_F; Fb_M];

    %% Restoring Forces
    % Described in Body Frame, at Center of Gravity. ALREADY FOR THE LEFT HAND
    % SIDE
    Fr_o = -[(W - B)*sin(theta);
                  -(W - B)*cos(theta)*sin(phi);
                  -(W - B)*cos(theta)*cos(phi);
                  -(rg(2)*W - rb(2)*B)*cos(theta)*cos(phi) + (rg(3)*W - rb(3)*B)*cos(theta)*sin(phi);
                   (rg(3)*W - rb(3)*B)*sin(theta) + (rg(1)*W - rb(1)*B)*cos(theta)*cos(phi);
                  -(rg(1)*W - rb(1)*B)*cos(theta)*sin(phi) - (rg(2)*W - rb(2)*B)*sin(theta)];

    %% Coriolis Forces
    % Described in Body Frame (specifically for sphere). ALREADY FOR THE LEFT
    % HAND SDE
    Crb = [0 -m*r m*q m*zg*r 0 0;
                m*r 0 -m*p 0 m*zg*r 0;
                -m*q m*p 0 -m*zg*p -m*zg*q 0;
                -m*zg*r 0 m*zg*p 0 I*r -I*q;
                0 -m*zg*r m*zg*q -I*r 0 I*p;
                0 0 0 I*q -I*p 0];  % Rigid Body Coriolis Force Matrix, at Center of Gravity
    Crb_o = Transform(Crb, rg_o);

    Ca = -ma * [0 0 0 0 w -v;
                  0 0 0 -w 0 u;
                  0 0 0 v -u 0;
                  0 w -v 0 0 0;
                  -w 0 u 0 0 0;
                  v -u 0 0 0 0];        % Added mass Coriolis Force Matrix, at Center of Buoyancy
    Ca_o = Transform(Ca, rb_o);

    Fc_o = (Crb_o + Ca_o) * V_t;            % Total Coriolis Force Matrix
    % Fc_0 = zeros(6,1);

    %% Damping Forces
    % Described in Body frame, at Center of Buoyancy
    % NOTE:(rg(3)*W - rb(3)*B)*sin(theta)
    % Normally, when we do velocity reading using sensor, it is the velocity at
    % the sensor location. Velocity has to be transformed first to the center
    % of buoyancy.
    % 
    % V_buoyancy = H(r_buoyancy_sensor) * V_sensor
    %
    % Compute the drag force, then transform the drag force matrix from
    % buoyancy center to the origin frame. FOR THE LEFT HAND SIDE

    % Drag Coefficient for sphere
    Re = abs(V_t(1:3)) * rho_fluid * D / mu_fluid;

    Cd11 = sphere_CD(Re(1));
    Cd22 = sphere_CD(Re(2));
    Cd33 = sphere_CD(Re(3));

    % Projected Surface Area
    Asp = pi * D^2;        % Stands for "Projected Area of the Surface"

    % Coefficient
    K11 = 1/2 * rho_fluid * Asp * Cd11;
    K22 = 1/2 * rho_fluid * Asp * Cd22;
    K33 = 1/2 * rho_fluid * Asp * Cd33;
    K44 = 1/64 * rho_fluid * D^5 * Cd11;
    K55 = 1/64 * rho_fluid * D^5 * Cd22;
    K66 = 1/64 * rho_fluid * D^5 * Cd33;

    Kd = -diag([K11 K22 K33 K44 K55 K66]);

    Fd_o = Transform(Kd, rb_o) * abs(V_t) .* V_t;

    %% External Forces (Weight + Buoyancy)
    Input_F = [0; 0; 0];           % Forces in x-axis, y-axis, and z-axis (Body Frame)
    F_Coord = [0; 0; R];      % External forces exerted on the top of the sphere, in line with the center of gravity

    Ex_Force = Command_Force(Input_F, F_Coord);
    impulse_time = 0.001;

    Ft_o = Ex_Force;

    %% Acceleration Computation 
    Acc_G = (MT) \ (Ft_o + (Fc_o + Fd_o + Fr_o));
    acc_z = Acc_G(3);
end