function Fc = sphere_Fc(V_t)
global Param
u = V_t(1);
v = V_t(2);
w = V_t(3);
p = V_t(4);
q = V_t(5);
r = V_t(6);

Param.Crb = [0 -Param.m*r Param.m*q Param.m*Param.zg*r 0 0;
       Param.m*r 0 -Param.m*p 0 Param.m*Param.zg*r 0;
       -Param.m*q Param.m*p 0 -Param.m*Param.zg*p -Param.m*Param.zg*q 0;
       -Param.m*Param.zg*r 0 Param.m*Param.zg*p 0 Param.I*r -Param.I*q;
       0 -Param.m*Param.zg*r Param.m*Param.zg*q -Param.I*r 0 Param.I*p;
       0 0 0 Param.I*q -Param.I*p 0];  % Rigid Body Coriolis Force Matrix, at Center of Gravity
Param.Crb_o = Transform(Param.Crb, Param.rg_o);

Param.Ca = -Param.ma * [0 0 0 0 w -v;
                  0 0 0 -w 0 u;
                  0 0 0 v -u 0;
                  0 w -v 0 0 0;
                  -w 0 u 0 0 0;
                  v -u 0 0 0 0];        % Added mass Coriolis Force Matrix, at Center of Buoyancy
Param.Ca_o = Transform(Param.Ca, Param.rb_o);

Fc = (Param.Crb_o + Param.Ca_o) * V_t;            % Total Coriolis Force Matrix
end