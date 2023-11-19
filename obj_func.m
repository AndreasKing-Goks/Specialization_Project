%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% optim                                                                   %
%                                                                         %              
% Initialize and Optimization                                             %
%                                                                         %
% Created:      18.11.2023	Andreas Sitorus                               %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [acc_z] = obj_func(X)
    %% Add Path
    % Current dir
    currentDir = fileparts(mfilename('fullpath'));

    % Add the 'Util' path
    utilPath = fullfile(currentDir, 'Util');
    addpath(utilPath);

    % Add the 'Optimization' path
    optimPath = fullfile(currentDir, 'Optimization');
    addpath(optimPath);

    %% Initialize Workspace
    percent_foam = X(1);
    percent_metal = X(2);
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
    acc = Sphere_Model(Ex_Force, Pos_N, Velo_B);
    acc_z = acc(3);
end