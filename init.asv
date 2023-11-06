%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init                                                                    %
%                                                                         %              
% Initialize workspace                                                    %
%                                                                         %
% Created:      27.09.2023	Andreas Sitorus                               %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
%% Initialize Workspace
percent_foam = 0.6;
percent_metal = 0.15;
Param = Sphere_Parameters(percent_foam, percent_metal);

%% Input Force in Body Frame
Input_F = [0; 0; 0];           % Forces in x-axis, y-axis, and z-axis (Body Frame)
F_Coord = [0; 0; Param.R];      % External forces exerted on the top of the sphere, in line with the center of gravity

Ex_Force = Command_Force(Input_F, F_Coord);
impulse_time = 0.001;

%% Initial Condition
Pos_N = Param.IC_Pos;
Velo_B = Param.IC_Velo;

%% Sphere Dynamic Model [FOR CHECKING]
Acc_G = Sphere_Model(Ex_Force, Pos_N, Velo_B)

%% HELP READING Acceleration result
% Forces defined in NED at first. then transformed to the body coordinate
% Thus, positive sign means downwards
% Postive acceleration means Negatively Buoyant
% Negative acceleration means Positively Buoyant

%% Check Status
% f_h = Param.foam_height
% m_h = Param.metal_height
% c_t_height = Param.c_t_height
% c_b_height = Param.c_b_height
% 
% V_foam = Param.V_foam
% V_metal = Param.V_metal
% V_cavity_t = Param.V_c_t
% V_cavity_b = Param.V_c_b
% V_cavity = Param.V_cavity
% V_cavity_TRUE = Param.V-(Param.V_foam + Param.V_metal)
% V_total_c =(Param.V_foam + Param.V_metal + Param.V_cavity)
% V_total = Param.V

% R = Param.R
% 
% c_f = Param.centroid_foam
% c_m = Param.centroid_metal
% V_HS = Param.V_HS
% C_HS = Param.centroid_HS
% 
% c_c_t = Param.centroid_cavity_t
% c_c_b = Param.centroid_cavity_b
% c_c = Param.centroid_cavity