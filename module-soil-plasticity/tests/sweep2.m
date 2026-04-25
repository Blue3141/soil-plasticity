%% sweep2.m — Bekker contact validation against analytical predictions
%
%  Sweeps impact velocity v0, runs MBDyn with module-soil-plasticity,
%  and compares simulated z_max, rebound velocity and peak force against
%  closed-form Bekker-Wong predictions.
%
%  Analytical reference (1-DOF rigid mass, 3 identical pads):
%    KBek   = kc/b + kphi
%    z_max  = [(n+1)*M*v0^2 / (2*A_tot*KBek)]^(1/(n+1))
%    e_rest = sqrt((n+1)*(1-kp)/2)
%    v_reb  = e_rest * v0
%    FN_max = A_tot * KBek * z_max^n
%
%  Sources:
%    Bekker (1969) "Introduction to Terrain-Vehicle Systems"
%    Wong (2008) "Theory of Ground Vehicles", 4th ed.
%    Jiang et al. (2021) Acta Astronautica — DEM validation of
%      Bekker-type models for lunar regolith (DOI:10.1016/j.actaastro.2021.04.030)
%    Li et al. (2013) J. Terramech. — granular contact, DEM vs Bekker
%      (DOI:10.1016/j.jterra.2013.03.001)

clear; clc; close all;

%% =========================================================================
%% 0. BASE SIMULATION PARAMETERS
%% =========================================================================
Simparam.t0 = 0.;
Simparam.tf = 1.5;
Simparam.dt = 5e-5;          % 0.05 ms — sufficient for stiff Bekker spring

Simparam.max_iterations     = 100;
Simparam.out_freq            = 200;    % output every 200 steps = 0.01 s
Simparam.simtol              = 1.e-6;
Simparam.derivatives_coeff   = 1.e-12;
Simparam.mass_multiplier     = 1;

Simparam.M_payload = 90.0;      % capsule dry mass [kg]
Simparam.H_payload = 0.900;     % CoM height above attachment ring [m]
Simparam.R_payload = 0.620;     % CoM radial offset [m]

Simparam.g = 1.62;              % lunar gravity [m/s^2]
Simparam.p0 = [0.0, 0.0, 0.01];% initial position (slight clearance above ground)
Simparam.v0 = [0.0, 0.0, 0.0]; % will be set per run in the sweep
Simparam.w0 = [0.0, 0.0, 0.0];

%% Bekker-Wong soil parameters  (loose lunar regolith — Bekker 1969, Table 7-3)
Simparam.kc        = 1400.;    % [N/m^(n+1)]
Simparam.kphi      = 820000.;  % [N/m^(n+2)]
Simparam.n_bek     = 1.0;      % pressure-sinkage exponent [-]
Simparam.b_pad     = 0.05;     % pad contact radius [m]  (b = smaller dim)
Simparam.k_plastic = 0.5;      % plastic ratio  z0 = kp*zmax  [-]
Simparam.h_hydro   = 0.;       % OFF — hydrodynamic disabled for clean Bekker validation
Simparam.VN0       = 0.3;      % reference velocity for tanh ramp [m/s]
Simparam.ZN0       = 0.02;     % spatial ramp thickness [m]
Simparam.mu_ground = 0.4;      % Coulomb friction coefficient [-]

%% Absorber
Simparam.KDF          = 7e-1;
Simparam.csi          = 0.2;
Simparam.absorber_type = 3;    % viscous + friction

Default_simparam = Simparam;

%% =========================================================================
%% 1. MATERIAL PROPERTIES
%% =========================================================================
% Type 1: main legs (CFRP, D40/d36)
k = 1; D = 40e-3; d = 36e-3;
props(k).E = 91e9; props(k).G = 35e9; props(k).rho = 1600;
props(k).A = pi/4*(D^2-d^2);
props(k).I = pi/64*(D^4-d^4); props(k).J = 2*props(k).I;

% Type 2: cross-braces (CFRP, D16/d14)
k = 2; D = 16e-3; d = 14e-3;
props(k).E = 91e9; props(k).G = 35e9; props(k).rho = 1600;
props(k).A = pi/4*(D^2-d^2);
props(k).I = pi/64*(D^4-d^4); props(k).J = 2*props(k).I;

% Type 3: main cross-members (CFRP, D27/d24)
k = 3; D = 27e-3; d = 24e-3;
props(k).E = 91e9; props(k).G = 35e9; props(k).rho = 1600;
props(k).A = pi/4*(D^2-d^2);
props(k).I = pi/64*(D^4-d^4); props(k).J = 2*props(k).I;

% Type 4: foot (Al, D47/d43)
k = 4; D = 47e-3; d = 43e-3;
props(k).E = 70e9; props(k).G = 26e9; props(k).rho = 2700;
props(k).A = pi/4*(D^2-d^2);
props(k).I = pi/64*(D^4-d^4); props(k).J = 2*props(k).I;

% Type 5: central column (Al, D47/d43)
k = 5; D = 47e-3; d = 43e-3;
props(k).E = 70e9; props(k).G = 26e9; props(k).rho = 2700;
props(k).A = pi/4*(D^2-d^2);
props(k).I = pi/64*(D^4-d^4); props(k).J = 2*props(k).I;

%% =========================================================================
%% 2. GEOMETRY  [mm, converted below]
%% =========================================================================
aste = [
    % main leg
    0, 0, 330,      1435, 0, 330,   1;
    % lower horizontal
    0, 0, 155,      625,  0, 155,   3;
    % vertical braces
    230, 0, 155,    230,  0, 330,   2;
    550, 0, 155,    550,  0, 330,   2;
    % diagonal braces
    0,   0, 155,    230,  0, 330,   2;
    405, 0, 155,    230,  0, 330,   2;
    405, 0, 155,    550,  0, 330,   2;
    1435,0, 225,    1195, 0, 330,   2;
    % transverse brace
    625,-100,155,   625, 100, 155,  2;
    % diagonal members to top
    625,-100,155,   1005, 0, 330,   3;
    625, 100,155,   1005, 0, 330,   3;
    % foot
    1435, 0, 0,     1435, 0, 330,   4;
    ];

l = size(aste, 1);

% Replicate 120° and 240° rotated legs
c = -0.5; s = 0.5*sqrt(3);
aste = [aste;
    c*aste(:,1)-s*aste(:,2), s*aste(:,1)+c*aste(:,2), aste(:,3), ...
    c*aste(:,4)-s*aste(:,5), s*aste(:,4)+c*aste(:,5), aste(:,6:end);
    c*aste(:,1)+s*aste(:,2),-s*aste(:,1)+c*aste(:,2), aste(:,3), ...
    c*aste(:,4)+s*aste(:,5),-s*aste(:,4)+c*aste(:,5), aste(:,6:end)];

% Add central column and inter-leg braces
aste = [aste;
    0, 0, 155,           0, 0, 610,           5;
    aste(3,4:6),         aste(3+l,4:6),       2;
    aste(3+l,4:6),       aste(3+2*l,4:6),     2;
    aste(3,4:6),         aste(3+2*l,4:6),     2];

% Convert mm → m; append rod length column
aste = [aste(:,1:6)*1e-3, aste(:,7)];
l = size(aste, 1);
aste = [aste, zeros(l,1)];
for i = 1:l
    aste(i,end) = norm(aste(i,1:3) - aste(i,4:6));
end

aste = split_mbdyn_rods(aste);

%% =========================================================================
%% 3. TOTAL MASS (for analytical formula)
%% =========================================================================
M_legs = 0;
for i = 1:size(aste,1)
    pid = aste(i,7);
    M_legs = M_legs + props(pid).rho * props(pid).A * aste(i,end);
end
M_total = Default_simparam.M_payload + M_legs;
fprintf('M_payload = %.1f kg,  M_legs = %.1f kg,  M_total = %.1f kg\n', ...
    Default_simparam.M_payload, M_legs, M_total);

%% =========================================================================
%% 4. ANALYTICAL BEKKER PREDICTIONS
%% =========================================================================
v0_values = [0.5, 1.0, 1.5, 2.0, 2.5];   % impact velocities [m/s]
N_v = numel(v0_values);

kc     = Default_simparam.kc;
kphi   = Default_simparam.kphi;
n      = Default_simparam.n_bek;
b_pad  = Default_simparam.b_pad;
kp     = Default_simparam.k_plastic;
A_pad  = pi * b_pad^2;         % area of one circular pad [m^2]
N_pads = 3;                    % three feet
A_tot  = N_pads * A_pad;
KBek   = kc/b_pad + kphi;      % combined Bekker stiffness [N/m^(n+2)]

% Maximum sinkage: energy balance  0.5*M*v0^2 = A_tot*KBek*z^(n+1)/(n+1)
z_max_anal = ((n+1) * 0.5 * M_total .* v0_values.^2 ./ (A_tot * KBek)) .^ (1/(n+1));

% Restitution coefficient (Wong 2008, §7.4)
e_rest     = sqrt((n+1)*(1-kp)/2);
v_reb_anal = e_rest * v0_values;

% Peak Bekker force at z_max
FN_peak_anal = A_tot * KBek .* z_max_anal .^ n;

fprintf('\n--- Analytical Bekker predictions ---\n');
fprintf('e_rest = %.4f\n', e_rest);
fprintf('%8s %10s %12s %12s\n', 'v0 [m/s]', 'z_max [m]', 'v_reb [m/s]', 'FN_max [N]');
for k = 1:N_v
    fprintf('%8.2f %10.4f %12.4f %12.1f\n', ...
        v0_values(k), z_max_anal(k), v_reb_anal(k), FN_peak_anal(k));
end

%% =========================================================================
%% 5. SWEEP — run one MBDyn simulation per v0
%% =========================================================================
mbdyn_exec = '/usr/local/mbdyn/bin/mbdyn';
base_filename = 'sim_sweep_bekker';

node_id_payload = 5000;   % capsule CoM node label

z_max_sim   = NaN(1, N_v);
v_reb_sim   = NaN(1, N_v);
FN_peak_sim = NaN(1, N_v);

for k = 1:N_v
    v0  = v0_values(k);
    tag = sprintf('v%04.0f', v0*1000);
    mbd_file = sprintf('%s_%s.mbd', base_filename, tag);
    mov_file = sprintf('%s_%s.mov', base_filename, tag);
    usr_file = sprintf('%s_%s.usr', base_filename, tag);

    Simparam_loc      = Default_simparam;
    Simparam_loc.v0   = [0., 0., -v0];
    Simparam_loc.p0   = [0., 0., 0.01];

    fprintf('\n[%d/%d] v0 = %.2f m/s  →  generating %s\n', k, N_v, v0, mbd_file);
    generate_mbdyn_truss_contact6(aste, props, Simparam_loc, mbd_file, 0);

    [status, cmdout] = system(sprintf( ...
        'LD_LIBRARY_PATH=/lib/x86_64-linux-gnu:/usr/lib/x86_64-linux-gnu %s %s', ...
        mbdyn_exec, mbd_file));

    if status ~= 0
        fprintf('  MBDyn FAILED:\n%s\n', cmdout);
        continue
    end

    %% -- Read .mov kinematics --
    try
        data_mov = readmatrix(mov_file, 'FileType', 'text');
    catch ME
        fprintf('  Cannot read .mov: %s\n', ME.message);
        continue
    end

    payload = data_mov(data_mov(:,1) == node_id_payload, :);
    if size(payload, 1) < 5
        fprintf('  payload node (%d) not found or too few rows in .mov\n', node_id_payload);
        continue
    end

    % .mov columns: label  x y z  phi1 phi2 phi3  vx vy vz  wx wy wz  ax ay az
    z_pay  = payload(:, 4);    % vertical position of capsule CoM
    vz_pay = payload(:, 10);   % vertical velocity

    %% -- Read .usr contact forces --
    foot_node_labels = [];
    FN_total_ts = [];   % total contact force per timestep

    if isfile(usr_file)
        try
            raw_usr = readmatrix(usr_file, 'FileType', 'text');
            rows_e  = raw_usr(:,1) == 9001;
            if any(rows_e)
                foot_node_labels = unique(raw_usr(rows_e, 2));
                n_feet = numel(foot_node_labels);
                FN_col = raw_usr(rows_e, 8);   % col 8 = FN
                Nt = floor(length(FN_col) / n_feet);
                if Nt > 0
                    FN_mat   = reshape(FN_col(1:n_feet*Nt), n_feet, Nt)';
                    FN_total_ts = sum(FN_mat, 2);
                    FN_peak_sim(k) = max(FN_total_ts);
                end
            end
        catch ME
            fprintf('  Warning reading .usr: %s\n', ME.message);
        end
    else
        fprintf('  .usr not found — force data skipped\n');
    end

    %% -- Sinkage: minimum z of any foot node (from .mov) --
    if ~isempty(foot_node_labels)
        z_foot_min = Inf;
        for fi = 1:numel(foot_node_labels)
            rows_f = data_mov(:,1) == foot_node_labels(fi);
            if any(rows_f)
                z_foot_min = min(z_foot_min, min(data_mov(rows_f, 4)));
            end
        end
        if isfinite(z_foot_min)
            z_max_sim(k) = -z_foot_min;   % positive sinkage
        end
    end

    % Fallback: infer sinkage from CoM depression if foot nodes unavailable
    if isnan(z_max_sim(k))
        z_max_sim(k) = -(min(z_pay) - z_pay(1));
    end

    %% -- Rebound velocity: first peak of vz_pay after maximum sinkage --
    [~, idx_zmin] = min(z_pay);
    if idx_zmin < numel(vz_pay)
        vz_after = vz_pay(idx_zmin:end);
        v_reb_sim(k) = max(vz_after);
    else
        v_reb_sim(k) = vz_pay(end);
    end

    fprintf('  z_max = %.4f m  (anal %.4f m)\n',  z_max_sim(k),   z_max_anal(k));
    fprintf('  v_reb = %.4f m/s (anal %.4f m/s)\n', v_reb_sim(k), v_reb_anal(k));
    if ~isnan(FN_peak_sim(k))
        fprintf('  FN_pk = %.1f N  (anal %.1f N)\n', FN_peak_sim(k), FN_peak_anal(k));
    end

    %% -- Cleanup --
    if isfile(mbd_file); delete(mbd_file); end
    if isfile(mov_file); delete(mov_file); end
    if isfile(usr_file); delete(usr_file); end
end

%% =========================================================================
%% 6. SAVE
%% =========================================================================
save('sweep_bekker_results.mat', ...
    'v0_values', 'z_max_sim', 'v_reb_sim', 'FN_peak_sim', ...
    'z_max_anal', 'v_reb_anal', 'FN_peak_anal', ...
    'e_rest', 'KBek', 'A_tot', 'M_total', 'Default_simparam');

%% =========================================================================
%% 7. PLOTS
%% =========================================================================

%% -- z_max vs v0 --
figure('Name', 'Maximum sinkage vs impact velocity');
plot(v0_values, z_max_anal*100, 'b-o', 'LineWidth', 2, 'DisplayName', ...
    'Bekker analytical  [Wong 2008]');
hold on; grid on;
plot(v0_values, z_max_sim*100, 'r--s', 'LineWidth', 2, 'DisplayName', ...
    'MBDyn simulation');
xlabel('Impact velocity  v_0  [m/s]');
ylabel('Maximum sinkage  z_{max}  [cm]');
title('Sinkage vs impact velocity — Bekker-Wong validation');
legend('Location', 'northwest');
text(v0_values(end)*0.7, max(z_max_anal)*100*0.15, ...
    sprintf('K_{Bek}=%.0f N/m^{n+2}  n=%.1f  k_p=%.2f', KBek, n, kp), ...
    'FontSize', 9, 'Color', 'k');

%% -- Rebound velocity vs v0 --
figure('Name', 'Rebound velocity vs impact velocity');
plot(v0_values, v_reb_anal, 'b-o', 'LineWidth', 2, 'DisplayName', ...
    sprintf('Analytical  (e=%.3f)', e_rest));
hold on; grid on;
plot(v0_values, v_reb_sim, 'r--s', 'LineWidth', 2, 'DisplayName', ...
    'MBDyn simulation');
% Reference line for perfectly elastic and fully plastic limits
plot(v0_values, v0_values, 'k:', 'DisplayName', 'elastic (e=1)');
plot(v0_values, zeros(size(v0_values)), 'k--', 'DisplayName', 'plastic (e=0)');
xlabel('Impact velocity  v_0  [m/s]');
ylabel('Rebound velocity  v_{reb}  [m/s]');
title('Rebound velocity — Bekker plasticity validation');
legend('Location', 'northwest');

%% -- Peak normal force vs v0 --
figure('Name', 'Peak normal force vs impact velocity');
plot(v0_values, FN_peak_anal/1e3, 'b-o', 'LineWidth', 2, 'DisplayName', ...
    'Analytical  F_{max}=A_{tot}K_{Bek}z_{max}^n');
hold on; grid on;
plot(v0_values, FN_peak_sim/1e3, 'r--s', 'LineWidth', 2, 'DisplayName', ...
    'MBDyn simulation (sum 3 pads)');
xlabel('Impact velocity  v_0  [m/s]');
ylabel('Peak normal force  F_{N,max}  [kN]');
title('Peak contact force — Bekker-Wong validation');
legend('Location', 'northwest');

%% -- Relative error plots --
figure('Name', 'Relative errors');
valid = ~isnan(z_max_sim) & z_max_sim > 0;

subplot(3,1,1); hold on; grid on;
if any(valid)
    err_z = (z_max_sim(valid) - z_max_anal(valid)) ./ z_max_anal(valid) * 100;
    bar(v0_values(valid), err_z);
end
xlabel('v_0 [m/s]'); ylabel('Error [%]');
title('Relative error: maximum sinkage');
yline(0, 'k--');

subplot(3,1,2); hold on; grid on;
valid_v = ~isnan(v_reb_sim) & v_reb_anal > 0;
if any(valid_v)
    err_v = (v_reb_sim(valid_v) - v_reb_anal(valid_v)) ./ v_reb_anal(valid_v) * 100;
    bar(v0_values(valid_v), err_v);
end
xlabel('v_0 [m/s]'); ylabel('Error [%]');
title('Relative error: rebound velocity');
yline(0, 'k--');

subplot(3,1,3); hold on; grid on;
valid_f = ~isnan(FN_peak_sim) & FN_peak_anal > 0;
if any(valid_f)
    err_f = (FN_peak_sim(valid_f) - FN_peak_anal(valid_f)) ./ FN_peak_anal(valid_f) * 100;
    bar(v0_values(valid_f), err_f);
end
xlabel('v_0 [m/s]'); ylabel('Error [%]');
title('Relative error: peak normal force');
yline(0, 'k--');

%% =========================================================================
%% AUXILIARY FUNCTIONS
%% =========================================================================
function new_aste = split_mbdyn_rods(aste)
all_points = [aste(:,1:3); aste(:,4:6)];
tol = 1e-6;
unique_nodes = unique(round(all_points/tol)*tol, 'rows', 'stable');
new_aste = [];
for i = 1:size(aste,1)
    p1 = aste(i,1:3); p2 = aste(i,4:6);
    prop_id = aste(i,7:end);
    L_orig  = norm(p2-p1);
    nodes_on_beam = [];
    for j = 1:size(unique_nodes,1)
        pn = unique_nodes(j,:);
        d1 = norm(pn-p1); d2 = norm(pn-p2);
        if d1 > tol && d2 > tol && abs((d1+d2)-L_orig) < tol
            nodes_on_beam = [nodes_on_beam; pn];
        end
    end
    if isempty(nodes_on_beam)
        new_aste = [new_aste; aste(i,:)];
        new_aste(end,end) = new_aste(end,end) * 0.1;
    else
        dists = sqrt(sum((nodes_on_beam-p1).^2, 2));
        [~, idx] = sort(dists);
        ordered = [p1; nodes_on_beam(idx,:); p2];
        for kk = 1:size(ordered,1)-1
            new_aste = [new_aste; ordered(kk,:) ordered(kk+1,:) prop_id];
        end
    end
end
end
