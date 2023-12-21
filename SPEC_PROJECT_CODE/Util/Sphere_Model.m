function [Acc_G] = Sphere_Model(Ex_Force, Pos_N, Velo_B)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sphere_Model()                                                          %
%                                                                         %
% Compute the dynamics of the Sphere given the external forces, position, %
% and velocity described at the body frame                                %
%                                                                         %
% Argument:                                                               %
% Ex_Force  : External forces directed at sphere, described on Body Frame %
% Pos_N     : Position described at NED frame. Dim(3x1).                  %
% Velo_B    : Velocity described at Body frame. Dim(3x1).                 %
%                                                                         %
% Output:                                                                 %
% Acc_G     : Total acceleration of the Sphere Model                      %
%                                                                         %
% Created:      27.09.2023	Andreas Sitorus                               %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Param
%% Position
% Described in NED frame
x = Pos_N(1,1);
y = Pos_N(2,1);
z = Pos_N(3,1);
phi = Pos_N(4,1);
theta = Pos_N(5,1);
psi = Pos_N(6,1);

% General Vector
T_ = [x; y; z];         % Translational Vector
R_ = [phi, theta, psi]; % Rotational Vector

%% Velocity
% Described in Origin frame
u = Velo_B(1,1);
v = Velo_B(2,1);
w = Velo_B(3,1);
p = Velo_B(4,1);
q = Velo_B(5,1);
r = Velo_B(6,1);

% General Vector
V_ = [u; v; w];     % Translation Velocity vector, read "Nu"
W_ = [p; q; r];     % Angular velocity vector, read "Omega"

V_t = [V_ ; W_];    % 6DOF Velocity

% %% Gravity Force
% %Described in Origin Frame
% Fg_F = [Param.W * sin(theta);
%         -Param.W * cos(theta) * sin(phi);
%         -Param.W * cos(theta) * cos(phi)];  % Gravity Forces
% Fg_M = zeros(3,1);                          % Gravity Moment
% 
% Fg = [Fg_F ; Fg_M];
% 
% %% Buoyancy Force
% %Described in Body Frame (Depends on the submerged volume)
% Fb_F = [-Param.B * sin(theta);
%         Param.B * cos(theta) * sin(phi);
%         Param.B * cos(theta) * cos(phi)];   % Buoyancy Forces
% lever_arm = -Param.rb + Param.rg;            
% Fb_M = SS_(lever_arm) * Fb_F;               % Buoyancy Moment
% 
% Fb = [Fb_F; Fb_M];

%% Restoring Forces
% Described in Body Frame, at Center of Gravity. ALREADY FOR THE LEFT HAND
% SIDE
Param.Fr_o = -[(Param.W - Param.B)*sin(theta);
      -(Param.W - Param.B)*cos(theta)*sin(phi);
      -(Param.W - Param.B)*cos(theta)*cos(phi);
      -(Param.rg(2)*Param.W - Param.rb(2)*Param.B)*cos(theta)*cos(phi) + (Param.rg(3)*Param.W - Param.rb(3)*Param.B)*cos(theta)*sin(phi);
      (Param.rg(3)*Param.W - Param.rb(3)*Param.B)*sin(theta) + (Param.rg(1)*Param.W - Param.rb(1)*Param.B)*cos(theta)*cos(phi);
      -(Param.rg(1)*Param.W - Param.rb(1)*Param.B)*cos(theta)*sin(phi) - (Param.rg(2)*Param.W - Param.rb(2)*Param.B)*sin(theta)];

%% Coriolis Forces
% Described in Body Frame (specifically for sphere). ALREADY FOR THE LEFT
% HAND SDE
Param.Crb = [0 -Param.m*r Param.m*q Param.m*Param.zg*r 0 0;
       Param.m*r 0 -Param.m*p 0 Param.m*Param.zg*r 0;
       -Param.m*q Param.m*p 0 -Param.m*Param.zg*p -Param.m*Param.zg*q 0;
       -Param.m*Param.zg*r 0 Param.m*Param.zg*p 0 Param.I*r -Param.I*q;
       0 -Param.m*Param.zg*r Param.m*Param.zg*q -Param.I*r 0 Param.I*p;
       0 0 0 Param.I*q -Param.I*p 0];  % Rigid Body Coriolis Force Matrix, at Center of Gravity
Param.Crb_o = Transform(Param.Crb, Param.rg_o);

Param.Ca = -Param.ma * [0 0 0 0 w -v;
                  0 0 0 -w 0 u;
                  0 0 0 v -u 0;
                  0 w -v 0 0 0;
                  -w 0 u 0 0 0;
                  v -u 0 0 0 0];        % Added mass Coriolis Force Matrix, at Center of Buoyancy
Param.Ca_o = Transform(Param.Ca, Param.rb_o);

Param.Fc_o = (Param.Crb_o + Param.Ca_o) * V_t;            % Total Coriolis Force Matrix
% Param.Fc_0 = zeros(6,1);

%% Damping Forces
% Described in Body frame, at Center of Buoyancy
% NOTE:(Param.rg(3)*Param.W - Param.rb(3)*Param.B)*sin(theta)
% Normally, when we do velocity reading using sensor, it is the velocity at
% the sensor location. Velocity has to be transformed first to the center
% of buoyancy.
% 
% V_buoyancy = H(r_buoyancy_sensor) * V_sensor
%
% Compute the drag force, then transform the drag force matrix from
% buoyancy center to the origin frame. FOR THE LEFT HAND SIDE

% Drag Coefficient for sphere
Re = abs(V_t(1:3)) * Param.rho_fluid * Param.D / Param.mu_fluid;

Cd11 = sphere_CD(Re(1));
Cd22 = sphere_CD(Re(2));
Cd33 = sphere_CD(Re(3));

% Projected Surface Area
Asp = pi * Param.D^2;        % Stands for "Projected Area of the Surface"

% Coefficient
K11 = 1/2 * Param.rho_fluid * Asp * Cd11;
K22 = 1/2 * Param.rho_fluid * Asp * Cd22;
K33 = 1/2 * Param.rho_fluid * Asp * Cd33;
K44 = 1/64 * Param.rho_fluid * Param.D^5 * Cd11;
K55 = 1/64 * Param.rho_fluid * Param.D^5 * Cd22;
K66 = 1/64 * Param.rho_fluid * Param.D^5 * Cd33;

Kd = -diag([K11 K22 K33 K44 K55 K66]);

Param.Fd_o = Transform(Kd, Param.rb_o) * abs(V_t) .* V_t;

%% External Forces (Weight + Buoyancy)
Param.Ft_o = Ex_Force;

%% Acceleration Computation 
Acc_G = (Param.MT) \ (Param.Ft_o + (Param.Fc_o + Param.Fd_o + Param.Fr_o));
end