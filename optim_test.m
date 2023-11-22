clear
clc
%% SET THE NEEDED DIRECTORIES
% Current dir
currentDir = fileparts(mfilename('fullpath'));

% Add the 'Util' path
utilPath = fullfile(currentDir, 'Util');
addpath(utilPath);

% Add the 'Util' path
optimPath = fullfile(currentDir, 'Optimization');
addpath(optimPath);

%% INITIALIZATION
% Starting Point
pf = 0.25;
pm = 0.25;

X = [pf;pm]

% SGD Parameter
func = @obj_func;
dfunc = @complex_step;
X0 = X;
Eg = 1e-10;
Ea = 1e-6;
Er = 1e-2;
maxIter = 500;
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
dfunargs.func = @acc_func;
dfunargs.funargs = funargs;

%% OPTIMIZE
[X,f,n,n2,n3,hist] = SGD(func, dfunc, X0, Eg, Ea, Er, maxIter, alpha_1, alpha_max, mu1, mu2,rho, maxMinorIter, X_Tol, funargs, dfunargs);