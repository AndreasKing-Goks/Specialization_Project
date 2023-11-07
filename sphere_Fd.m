function Fd = sphere_Fd(V_t)
global Param
% Drag Coefficient for sphere
Re = V_t(1:3) * Param.rho_fluid * Param.D / Param.mu_fluid;

Cd11 = sphere_CD(Re(1));
Cd22 = sphere_CD(Re(2));
Cd33 = sphere_CD(Re(3));

% Projected Surface Area
Asp = pi * Param.D^2;        % Stands for "Projected Area of the Surface"

% Coefficient
K11 = 1/2 * Param.rho_fluid * Asp * Cd11;
K22 = 1/2 * Param.rho_fluid * Asp * Cd22;
K33 = 1/2 * Param.rho_fluid * Asp * Cd33;
K44 = 1/64 * Param.rho_fluid * Param.D^5 * Cd11;
K55 = 1/64 * Param.rho_fluid * Param.D^5 * Cd22;
K66 = 1/64 * Param.rho_fluid * Param.D^5 * Cd33;

Kd = -diag([K11 K22 K33 K44 K55 K66]);

Fd = Transform(Kd, Param.rb_o) * (V_t .* abs(V_t));
end