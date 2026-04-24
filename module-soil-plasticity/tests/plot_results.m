clc; clear; close all;

base_name = 'lunar_landing';

% =========================================================================
% Kinematics from .mov
% Columns: label  x y z  phi1 phi2 phi3  vx vy vz  wx wy wz  ax ay az  ...
% =========================================================================
mov_file = [base_name '.mov'];
data = readmatrix(mov_file, 'FileType', 'text');
node = data(data(:,1) == 1, :);

g    = -1.62;          % lunar gravity [m/s^2]
step = 1e-3;           % time step [s]

z  = node(:, 4);       % z position
vz = node(:, 10);      % z velocity
az = node(:, 16);      % z acceleration
t  = step * (1:length(z))';

% =========================================================================
% Contact data from .usr
% Columns: elem_label  node_label  Fx Fy Fz  F_BEK  FNH  FN  z0
% Written by SurfaceImpactElem::Output() -> OH.Loadable()
% =========================================================================
usr_file = [base_name '.usr'];

has_usr = isfile(usr_file);
if has_usr
    raw  = readmatrix(usr_file, 'FileType', 'text');
    rows = raw(:,1) == 1 & raw(:,2) == 1;
    if ~any(rows)
        warning('No rows matching elem=1 node=1 in %s', usr_file);
        has_usr = false;
    else
        F_BEK = raw(rows, 6);
        FNH   = raw(rows, 7);
        FN    = raw(rows, 8);
        z0    = raw(rows, 9);
        tp    = step * (1:length(F_BEK))';
    end
else
    warning('.usr file not found: %s  --  contact plots skipped', usr_file);
end

% =========================================================================
% Kinematics plots
% =========================================================================

figure('Name', 'Position');
plot(t, z); grid on;
xlabel('t [s]'); ylabel('z [m]');
title('Vertical position');

figure('Name', 'Velocity');
plot(t, vz); grid on;
xlabel('t [s]'); ylabel('v_z [m/s]');
title('Vertical velocity');

figure('Name', 'Net acceleration');
plot(t, az - g); grid on;
xlabel('t [s]'); ylabel('a_z - g  [m/s^2]');
title('Net acceleration (gravity removed)');

figure('Name', 'Acceleration vs sinkage');
plot(-z, az - g); grid on;
xlabel('penetration [m]'); ylabel('a_z - g  [m/s^2]');
title('Force-sinkage (via acceleration)');

figure('Name', 'Mechanical energy');
E = 0.5 * 500 * vz.^2 - 500 * g * (z - min(z));
plot(t, E); grid on;
xlabel('t [s]'); ylabel('E [J]');
title('Mechanical energy (kinetic + gravitational PE)');

% =========================================================================
% Contact force plots  (only if .usr was read successfully)
% =========================================================================
if ~has_usr
    return
end

figure('Name', 'Force breakdown');
hold on; grid on;
plot(tp, F_BEK, 'b',   'DisplayName', 'F_{BEK}  (Bekker)');
plot(tp, FNH,   'r',   'DisplayName', 'F_{NH}  (hydro)');
plot(tp, FN,    'k--', 'DisplayName', 'F_N  (total)');
xlabel('t [s]'); ylabel('F [N]');
title('Normal force components');
legend('Location', 'best');

figure('Name', 'Bekker force');
plot(tp, F_BEK); grid on;
xlabel('t [s]'); ylabel('F_{BEK} [N]');
title('Bekker-Wong normal force');

figure('Name', 'Hydrodynamic force');
plot(tp, FNH); grid on;
xlabel('t [s]'); ylabel('F_{NH} [N]');
title('Hydrodynamic normal force');

figure('Name', 'Total normal force');
plot(tp, FN); grid on;
xlabel('t [s]'); ylabel('F_N [N]');
title('Total normal force');

figure('Name', 'Permanent sinkage');
plot(tp, z0); grid on;
xlabel('t [s]'); ylabel('z_0 [m]');
title('Permanent sinkage (Bekker plasticity)');

% Hysteresis loop: Bekker force vs penetration
sinkage = interp1(t, -z, tp, 'linear', 0);
figure('Name', 'Bekker hysteresis');
plot(sinkage, F_BEK); grid on;
xlabel('penetration z [m]'); ylabel('F_{BEK} [N]');
title('Bekker load-unload hysteresis');

% Force ratio: hydro vs total
figure('Name', 'Hydro fraction');
frac = zeros(size(FN));
ok   = FN > 0;
frac(ok) = FNH(ok) ./ FN(ok) * 100;
plot(tp, frac); grid on;
xlabel('t [s]'); ylabel('F_{NH} / F_N  [%]');
title('Hydrodynamic fraction of total normal force');
