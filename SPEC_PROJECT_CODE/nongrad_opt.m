clear
clc
%% SET THE NEEDED DIRECTORIES
% Current dir
currentDir = fileparts(mfilename('fullpath'));

% Add the 'Util' path
utilPath = fullfile(currentDir, 'Util');
addpath(utilPath);

% Add the 'Util' path
optimPath = fullfile(currentDir, 'Non_Gradient_Based');
addpath(optimPath);

%% INITIALIZATION
% Starting Point
pf = 0.6;
pm = 0.6;

X = [pf;pm]';
% X = rand(1, 2);

% Nedler Mead Parameter
func = @obj_func;
X0 = X;
alpha = 1; 
beta = 2; 
gamma = 1/2; 
delta = 1/2;
maxIter = 50;
Ea = 1e-6;
Er = 1e-2;
X_Tol = 1e-20;
step_size = 0.05; 

%% OPTIMIZE
[traces, opt_point, opt_value] = nelder_mead(func, X0, alpha, beta, gamma, delta, maxIter, Ea, Er, X_Tol, step_size);
fprintf('Total system acceleration at optimum point is:%f m/s2\n',acc_func(opt_point))

%% STEP PLOT

% Create meshgrid
% max_pf = 0.9;
% min_pf = 0.2;
% 
% max_pm = 0.5;
% min_pm = 0.25;

max_pf = 1.0;
min_pf = 0.0;

max_pm = 1.0;
min_pm = 0.0;

resolution = 0.01;

percent_foam = min_pf:resolution:max_pf; 
percent_metal = min_pm:resolution:max_pm;
size_mesh = [length(percent_foam), length(percent_metal)];
[pf, pm] = meshgrid(percent_foam, percent_metal);

% Plot result
% Initialize an array to store the z-acceleration values
acc_z = zeros(size_mesh);
X = zeros(2,1);

% Calculate acceleration for each point in the mesh grid
for row = 1:size(acc_z, 1)
    for column = 1:size(acc_z, 2)
        % Initialize parameters
        Param = Sphere_Parameters(percent_foam(row), percent_metal(column));

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
figure('Name','SGD Plot','NumberTitle','off');
contour(pf, pm, acc_z', 200);
hold on

% Plot the descent point
scatter(traces(:,1), traces(:,2), 'r', 'filled');
hold on

% Plot the nelder triangle
for i = 1:size(traces, 1)-2
    vertices = traces(i:i+2, 1:2); % Extract the x and y coordinates of the vertices of the triangle
    point1 = vertices(1, :);
    point2 = vertices(2, :);
    point3 = vertices(3, :);
    
    plot([point1(1), point2(1)], [point1(2), point2(2)], 'k-'); % Line from point1 to point2
    hold on;
    plot([point2(1), point3(1)], [point2(2), point3(2)], 'k-'); % Line from point2 to point3
    hold on;
    plot([point3(1), point1(1)], [point3(2), point1(2)], 'k-'); % Line from point3 to point1
    hold on;
end

xlabel(['Floater height (as a percentage of radius R = ' num2str(Param.R) ') meters']);
ylabel(['Weight height (as a percentage of radius R = ' num2str(Param.R) ') meters']);
title('Acceleration Contour Plot');
colorbar; % Add a colorbar to indicate the acceleration values

%% PLOT VALUE DECREASE
figure('Name','Value Decrease Plot','NumberTitle','off');
plot(traces(:,3));
title('Acceleration Value Decrease')
xlabel('Iteration')
ylabel('Acceleration value [m/s2]')
grid on