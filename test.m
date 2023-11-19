clear
clc

% Starting Point
pf = 0.47;
pm = 0.145;

X = [pf;pm];

% Declare function arguments
funargs = [];

% Declare derivative function arguments
dfunargs.obj_func = @obj_func;
dfunargs.funargs = funargs;

% Compute gradient
acc_z = obj_func(X)
grad_fd = finite_difference(X, dfunargs)