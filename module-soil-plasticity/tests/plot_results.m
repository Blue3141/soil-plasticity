clc; clear; close all;

base_name = 'lunar_landing';

mov_file = [base_name '.mov'];
log_file = [base_name '.log'];

% ---- kinematics from .mov ----
data = readmatrix(mov_file, 'FileType', 'text');
node = data(data(:,1)==1,:);

g    = -1.62;
step = 1e-3;

z  = node(:,4);
vz = node(:,10);
az = node(:,16);
t  = step * (1:length(z))';

% ---- contact data from .log ----
% Columns: elem_label  node_label  Fx Fy Fz  F_BEK  FNH  FN  z0
raw  = readmatrix(log_file, 'FileType', 'text');
rows = raw(:,1)==1 & raw(:,2)==1;
F_BEK = raw(rows, 6);
FNH   = raw(rows, 7);
FN    = raw(rows, 8);
z0    = raw(rows, 9);
tp    = step * (1:length(F_BEK))';

% ======== kinematics ========

figure('Name','Position');
plot(t, z); grid on;
xlabel('t [s]'); ylabel('z [m]');

figure('Name','Velocity');
plot(t, vz); grid on;
xlabel('t [s]'); ylabel('v_z [m/s]');

figure('Name','Net acceleration');
plot(t, az - g); grid on;
xlabel('t [s]'); ylabel('a_z - g [m/s^2]');

figure('Name','Acceleration vs sinkage');
plot(-z, az - g); grid on;
xlabel('penetration [m]'); ylabel('a_z - g [m/s^2]');

figure('Name','Mechanical energy');
E = 0.5*500*vz.^2 - 500*g*(z - min(z));
plot(t, E); grid on;
xlabel('t [s]'); ylabel('E [J]');

% ======== contact forces ========

figure('Name','Force breakdown');
hold on; grid on;
plot(tp, F_BEK, 'b',  'DisplayName','F_{BEK}');
plot(tp, FNH,   'r',  'DisplayName','F_{NH}');
plot(tp, FN,    'k--','DisplayName','F_N total');
xlabel('t [s]'); ylabel('F [N]');
legend('Location','best');

figure('Name','Bekker force');
plot(tp, F_BEK); grid on;
xlabel('t [s]'); ylabel('F_{BEK} [N]');

figure('Name','Hydrodynamic force');
plot(tp, FNH); grid on;
xlabel('t [s]'); ylabel('F_{NH} [N]');

figure('Name','Total normal force');
plot(tp, FN); grid on;
xlabel('t [s]'); ylabel('F_N [N]');

figure('Name','Permanent sinkage');
plot(tp, z0); grid on;
xlabel('t [s]'); ylabel('z_0 [m]');

% force vs sinkage (hysteresis loop)
sinkage = interp1(t, -z, tp, 'linear', 0);
figure('Name','Bekker hysteresis');
plot(sinkage, F_BEK); grid on;
xlabel('penetration z [m]'); ylabel('F_{BEK} [N]');
