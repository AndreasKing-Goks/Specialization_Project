%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% obj_func_plot                                                           %
%                                                                         %              
% Plot the objective function contour plot                                %
%                                                                         %
% Created:      18.11.2023	Andreas Sitorus                               %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
%% Add Path
% Current dir
currentDir = fileparts(mfilename('fullpath'));

% Add the 'Util' path
utilPath = fullfile(currentDir, 'Util');
addpath(utilPath);

%% Create meshgrid
percent_foam_r = 0.0:0.005:1.0; 
percent_metal_r = 0.0:0.005:1.0;
size_mesh = [length(percent_foam_r), length(percent_metal_r)];
[pf, pm] = meshgrid(percent_foam_r, percent_metal_r);

%% Plot result
% Initialize an array to store the z-acceleration values
acc_z = zeros(size_mesh);

% Calculate acceleration for each point in the mesh grid
for row = 1:size(acc_z, 1)
    for column = 1:size(acc_z, 2)
        % Initialize parameters
        Param = Sphere_Parameters(pf(row, column), pm(row, column));

        % Initial condition
        Pos_N = Param.IC_Pos;
        Velo_B = Param.IC_Velo;

        % Initial force in body-fixed frame
        Input_F = [0; 0; 0];           % Forces in x-axis, y-axis, and z-axis (Body Frame)
        F_Coord = [0; 0; Param.R];     % External forces exerted on the top of the sphere, in line with the center of gravity

        Ex_Force = Command_Force(Input_F, F_Coord);
        impulse_time = 0.001;

        % Compute Acceleration
        acc = Sphere_Model(Ex_Force, Pos_N, Velo_B);
        acc_z(row,column) = abs(acc(3));
    end
end

% Create a contour plot
contour(pf, pm, acc_z, 500); 
xlabel(['Floater height (as a percentage of radius R = ' num2str(Param.R) ' meters)']);
ylabel(['Weight height (as a percentage of radius R = ' num2str(Param.R) ') meters']);
title('Acceleration Contour Plot');
colorbar; % Add a colorbar to indicate the acceleration values

%% HELP READING Acceleration result
% Forces defined in NED at first, then transformed to the body coordinate
% Thus, positive sign means downwards
% Postive acceleration means Negatively Buoyant
% Negative acceleration means Positively Buoyant