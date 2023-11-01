clear
clc
%% PLOT
data = load('Result/check_acceleration.mat');
Acc = data.Acc_B_S;
Velo = data.Velo_B_S;
Pos = data.Pos_N_S;
tout = data.tout;

Re = zeros(length(tout),3);
Cd = zeros(length(tout),3);

Velo_T = Velo.Data(:,1:3);

%% Iterate over Velo_T
rho_fluid = 1000;
D = 0.25;
mu_fluid = 1.002e-3;

[n,m] = size(Velo_T);
for col = 1:m
    for row = 1:n
        Rn =  Velo_T(row,col) * rho_fluid * D / mu_fluid;
        Re(row,col) = Rn;
        Cd(row,col) = sphere_CD(Rn);
    end
end

%% PLOT
% figure(1)
% plot(tout(1:4001), Pos.Data(1:4001,3));               
% xlabel('Water depth (meter)');          
% ylabel('Time (second)');          
% title('Position evolution');
% grid on
% 
% figure(2)
% plot(tout(1:4001), Velo.Data(1:4001,3));               
% xlabel('Velocity (m/s)');          
% ylabel('Time (second)');          
% title('Velocity evolution');
% grid on

figure(3)
plot(tout(1:101), Acc.Data(1:101,3));               
xlabel('Acceleration (m/s2)');          
ylabel('Time (second)');          
title('Acceleration evolution');
grid on

figure(4)
plot(tout(1:101), Re(1:101,3));               
xlabel('Reynold Number');          
ylabel('Time (second)');          
title('Reynold Number evolution');
grid on

figure(5)
plot(tout(1:101), Cd(1:101,3));               
xlabel('Cd');          
ylabel('Time (second)');          
title('Drag Coefficient evolution');
grid on

