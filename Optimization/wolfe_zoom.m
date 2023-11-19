function [alph_s, n] = wolfe_zoom(alph_lo, alph_hi, fun, dfun, pk, x0, mu1,...
    mu2, maxIter,funargs,dfunargs)
% function to carry out the pinpoint (or zoom) part of the line search
% this is called by strongwolfe

% Some initializations:
done = 0;
fk = fun(x0,funargs);
gk = dfun(x0,dfunargs);
alph = zeros(1,maxIter);
falph = zeros(1,maxIter);
dfalph = zeros(1,maxIter);

n = 1;
flo = fun(x0+alph_lo*pk,funargs); % function evaluation at alpha low

while (~done && n<maxIter)
    alph(n) = (alph_lo + alph_hi)/2; % new intermediate alpha value by 
    %simple interpolation
    falph(n) = fun(x0+alph(n)*pk,funargs); % evaluate f a new point
    
    if (falph(n) > fk + mu1*alph(n)*(gk'*pk)) || (falph(n) > flo)
        alph_hi = alph(n);
    else
        dfalph(n) = dfun(x0+alph(n)*pk,dfunargs).'*pk;
        if abs(dfalph(n)) <= mu2*abs(gk'*pk)
            done = 1;
        elseif dfalph(n)*(alph_hi-alph_lo) >=0
            alph_hi = alph(n);
        else
            alph_lo = alph(n);
        end
    end
    n = n+1;
    
end

alph_s = alph(n-1);