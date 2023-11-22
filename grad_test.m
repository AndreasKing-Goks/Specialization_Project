clear
clc

% Current dir
currentDir = fileparts(mfilename('fullpath'));

% Add the 'Util' path
utilPath = fullfile(currentDir, 'Util');
addpath(utilPath);

% Add the 'Util' path
optimPath = fullfile(currentDir, 'Optimization');
addpath(optimPath);

% Starting Point
pf = 0.47;
pm = 0.145;

X = [pf;pm]

% Declare function arguments
funargs = [];

% Declare derivative function arguments
dfunargs.func = @acc_func;
dfunargs.funargs = funargs;

% Compute acceleration
acc_z = acc_func(X)

% Compute gradient
grad_cs = complex_step(X, dfunargs)
grad_fd = finite_difference(X, dfunargs)