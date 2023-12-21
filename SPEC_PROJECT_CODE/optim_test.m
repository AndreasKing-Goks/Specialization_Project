clear
clc
%% SET THE NEEDED DIRECTORIES
% Current dir
currentDir = fileparts(mfilename('fullpath'));

% Add the 'Util' path
utilPath = fullfile(currentDir, 'Util');
addpath(utilPath);

% Add the 'Util' path
optimPath = fullfile(currentDir, 'Gradient_Based');
addpath(optimPath);

%% INITIALIZATION
% Starting Point
pf = 0.1;
pm = 0.9;

% X = [pf;pm]
X = rand(2, 1);

% SGD Parameter
func = @acc_func;
dfunc = @complex_step;
X0 = X;
Eg = 1e-10;
Ea = 1e-6;
Er = 1e-2;
maxIter = 40;
alpha_1 = 0.01;
alpha_max = 0.05;
mu1 = 1e-4;
mu2 = 1e-6;
rho = 2;
maxMinorIter = 20;
X_Tol = 1e-6;

% Declare function arguments
funargs = [];

% Declare derivative function arguments
dfunargs.func = func;
dfunargs.funargs = funargs;

%% OPTIMIZE
[X,f,n,n2,n3,hist] = SGD(func, dfunc, X0, Eg, Ea, Er, maxIter, alpha_1, alpha_max, mu1, mu2,rho, maxMinorIter, X_Tol, funargs, dfunargs);

%% STEP PLOT
% Unpack history
x1 = hist(:,1);
x2 = hist(:,2);
z = hist(:,3);

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
        acc_z(row,column) = acc(3);
    end
end

% Create a contour plot
figure('Name','SGD Plot','NumberTitle','off');
contour(pf, pm, acc_z', 200);
hold on
plot(x1,x2, 'k')
hold on
scatter(x1,x2, 'k')
hold on

xlabel(['Floater height (as a percentage of radius R = ' num2str(Param.R) ') meters']);
ylabel(['Weight height (as a percentage of radius R = ' num2str(Param.R) ') meters']);
title('Acceleration Contour Plot');
colorbar; % Add a colorbar to indicate the acceleration values

%% PLOT VALUE DECREASE
figure('Name','Value Decrease Plot','NumberTitle','off');
plot(z);
title('Acceleration Value Decrease')
xlabel('Iteration')
ylabel('Acceleration value [m/s2]')
grid on

%% PLOT MINOR PHASE
figure('Name','Line Search Counter','NumberTitle','off');
stem(n2);
hold on
stem(n3);
hold on
title('Line Search Phase Counter Over Time ')
xlabel('Iteration')
ylabel('Count')
legend('Strong Wolfe Condition','Zooming Phase')
grid on
