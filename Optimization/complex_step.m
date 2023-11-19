function [grad] = complex_step(X, dfunargs)
    % Unpack derivative arugments
    obj_func = dfunargs.obj_func;
    funargs = dfunargs.funargs;

    % Step size
    h = 1e-10;

    % Number of X elements
    n = numel(X);

    % Gradient memory
    grad = zeros(size(X,1), size(X,2));

    % Copy of original X vector
    X_complex = X;
    
    % Looping process
    for i = 1:n
        % Do complex step
        X_complex(i) = X_complex(i) + 1i * h;
        
        % Compute complex step gradient
        if nargin(obj_func) == 1
            grad(i) = imag(obj_func(X_complex)) / h;
        else
            grad(i) = imag(obj_func(X_complex, funargs)) / h;
        end
    end
end