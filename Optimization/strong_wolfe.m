function [x,f,n,n2,alph_s] = strong_wolfe(fun, dfun, pk, x0, alph_1, alph_max, mu1, mu2, rho, maxIter,funargs,dfunargs)
n = 1; % start counting major line search iterations
done = 0; % flag for return
alph = zeros(1,maxIter); % initialize alpha vector
alph(2) = alph_1; % initial step
phi = zeros(1,maxIter); % initialize phi vector
dphi = zeros(maxIter); % initialize dphi
n2 = 0; % counter for zoom stage

% find the function value at the starting point
phi0 = fun(x0,funargs); % note that the function must allow multiple inputs
% find the gradient at the starting point
phiPrime0 = dfun(x0,dfunargs).'*pk;  % note that the function must allow multiple inputs


while (~done && n<maxIter) % loop until we reach max iterations or we are happy

    phi(n) = fun(x0+alph(n)*pk,funargs); % evaluate the function 

    if ( (phi(n) > phi0 + mu1*alph(n)*(phiPrime0)) || ( n>1 && phi(n) > phi(n-1) ) ) 
        % Then the sufficient decrease is not fullfilled, but the minimum
        % is bracketed
        if n>1
            [alph_s, n2] = wolfezoom(alph(n-1), alph(n), fun, dfun, pk, x0, mu1, mu2, maxIter,funargs,dfunargs );
        else
            [alph_s, n2] = wolfezoom(0,alph_1, fun, dfun, pk, x0, mu1, mu2, maxIter,funargs,dfunargs);
        end

        x = x0+alph_s*pk;
        f = fun(x,funargs);
        done = 1;
    else % Then the sufficient decrease is fullfilled
        dphi(n) = (dfun(x0+alph(n)*pk,dfunargs)).'*pk;
        if abs(dphi(n)) <= -mu2*phiPrime0 % Then the curvature condition is fullfilled
            alph_s = alph(n);
            x = x0+alph_s*pk;
            f = fun(x,funargs);
            done = 1;    
        elseif dphi(n) >= 0
            if n>1
                [alph_s, n2] = wolfezoom(alph(n),alph(n-1), fun, dfun, pk, x0, mu1, mu2, maxIter,funargs,dfunargs);
            else
                [alph_s, n2] = wolfezoom(0,alph(n), fun, dfun, pk, x0, mu1, mu2, maxIter,funargs,dfunargs);
            end
            x = x0+alph_s*pk;
            f = fun(x, funargs);
            done = 1;
        else
            alph(n+1) = min([alph(n)*rho, alph_max]); % look farther
            if alph(n+1) >= alph_max
                alph_s = alph_max; 
                x = x0+alph_s*pk;
                f = fun(x,funargs);
                done = 1;
            end
        end
    end
    n = n+1; 
end

if n==maxIter
    % max Wolfe bracket iterations
    alph_s=alph(end);
    x = x0+alph_s*pk;
    f = fun(x,funargs); 
end