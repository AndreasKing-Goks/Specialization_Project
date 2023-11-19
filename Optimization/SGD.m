function [x,f,n,n2,n3,hist] = SGD(fun, dfun, x0, Eg, Ea, Er, maxIter, alph_1, ...
    alph_max, mu1, mu2,rho, maxMinorIter, xTol,funargs,dfunargs)

% Minimize f using the steepest descent method
% returns the vector x, minimum f, number of iterations n, strong Wolfe
% iterations n2, zoom iterations n3, and the search history (x,f) for each
% major iteration

% input: function fun, function dfun for the gradient, starting point x0,
% convergence parameters Ea, Er, Eg, the maximum number of major iterations maxIter,
% and strong Wolfe search parameters

n2 = zeros(maxIter,1);
n3 = zeros(maxIter,1);
hist = 0;
k = 1;
xs = zeros(length(x0),maxIter+1);
fs = zeros(maxIter+1,1);
cond1 = zeros(maxIter+1,1);
gk = zeros(length(x0),maxIter);
pk = zeros(length(x0),maxIter); 
done = 0;
xchange = xTol+1;
xs(:,1) = x0;

while k<=maxIter && ~done
    % compute the gradient of the function
    gk(:,k) = dfun(xs(:,k),dfunargs);
    
    % check for convergence
    if k>1
        xchange = abs(xs(:,k) - xs(:,k-1) );
    end
    

    if (norm(gk(:,k),2)<= Eg && ( max(xchange)<xTol))
            x = xs(:,k);
            f = fun(x,funargs);
            done = 1;
    else
        
        pk(:,k) = -gk(:,k)/norm(gk(:,k),2);

        % compute a line search to find the step alph_s in the direction of pk
        [xs(:,k+1),fs(k+1),n2(k),n3(k),alph_s] = strongwolfe(fun, dfun, pk(:,k), xs(:,k), alph_1,...
            alph_max, mu1, mu2,rho, maxMinorIter,funargs,dfunargs);
        % store alph_s as the step size for the next search
        if k>1
            alph_1 = alph_s*(gk(:,k-1).'*pk(:,k-1))/(gk(:,k).'*pk(:,k)); 
        else
            alph_1 = alph_s; 
        end

        fs(k) = fun(xs(:,k),funargs);

        cond1(k) = abs(fs(k+1)-fs(k)) <= (Ea+Er*abs(fs(k)));
        % check for convergence
        if (k>1) && (cond1(k) == 1) 
            xchange = abs(xs(:,k) - xs(:,k-1) );
            if max(xchange)<xTol
                x = xs(:,k);
                f = fs(k);
                done = 1;
            end
        end
    end
    
    k = k+1;
end

if k>=maxIter
    x = xs(:,k);
    f = fs(k);
end
n = k-1;
hist = [xs(:,1:n).' fs(1:n)];
n2 = n2(1:n); 
n3 = n3(1:n); 