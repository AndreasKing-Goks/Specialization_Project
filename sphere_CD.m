function CD = sphere_CD(Re)
% This function computes the drag coefficient of a sphere as a function of 
% the Reynolds number Re.
% Curve fitted after fig . A -56 in Evett and Liu: 
% "Fluid Mechanics and Hydraulics"
% Website: https://folk.ntnu.no/leifh/teaching/tkt4140/._main012.html
    if Re <= 0.0
        CD = 0.0;
    elseif Re > 8.0e6
        CD = 0.2;
    elseif Re > 0.0 && Re <= 0.5
        CD = 24.0/Re;
    elseif Re > 0.5 && Re <= 100.0
        p = [4.22 -14.05 34.87 0.658];
        CD = polyval(p, 1.0/Re);
    elseif Re > 100.0 && Re <= 1.04e4
        p = [-30.41 43.72 17.08 2.41];
        CD = polyval(p, 1.0/log10(Re));
    elseif Re > 1.04e4 && Re <= 3.35e5
        p = [-0.1584 2.031 -8.472 11.932];
        CD = polyval(p, log10(Re));
    elseif Re > 3.35e5 && Re <= 5.0e5
        x1 = log10(Re/4.5e5);
        CD = 91.08*x1^4 + 0.0764;
    else
        p = [-0.06338 1.1905 -7.332 14.93];
        CD = polyval(p, log10(Re));
    end
end