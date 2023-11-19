function [grad] = finite_difference(X, dfunargs)
    % Unpack derivative arugments
    obj_func = dfunargs.obj_func;
    funargs = dfunargs.funargs;

    % Step size
    h = 1e-5;

    % Number of X elements
    n = numel(X);

    % Gradient memory
    grad = zeros(size(X,1), size(X,2));

    % Copy of the original vector X
    X_perturbed = X;
    
    % Looping process
    for i = 1:n
        % Perturb the X vector
        X_perturbed(i) = X_perturbed(i) + h;
        
        % Compute complex step gradient
        if nargin(obj_func) == 1
            grad(i) = (obj_func(X_perturbed)-obj_func(X)) / h;
        else
            grad(i) = (obj_func(X_perturbed, funargs)-obj_func(X, funargs)) / h;
        end
    end
end