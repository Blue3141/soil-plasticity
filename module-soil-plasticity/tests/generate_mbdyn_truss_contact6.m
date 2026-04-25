function generate_mbdyn_truss_contact6(rod_data, properties, Simparam, filename, recover)
% generate_mbdyn_truss_contact6 — write an MBDyn input file for a 3-leg
% lander truss using the surface-impact Bekker-Wong soil module.
%
% Replaces the old cont-contact rod chain (shadow nodes 4000+, extra nodes
% 9000+, LuGre friction joints) with a single UserDefinedElem that acts
% directly on the three foot nodes.
%
% Soil parameters expected in Simparam:
%   .kc         Bekker cohesive modulus          [N/m^(n+1)]
%   .kphi       Bekker frictional modulus         [N/m^(n+2)]
%   .n_bek      Sinkage exponent                  [-]
%   .b_pad      Pad radius (smaller dimension)    [m]
%   .k_plastic  Plastic ratio k_p in (0,1)        [-]
%   .h_hydro    Hydrodynamic coeff h=0.5*rho*Cd   [kg/m^3]; 0 = disabled
%   .VN0        Velocity scale for tanh stab.     [m/s]
%   .ZN0        Depth ramp scale from z0          [m]
%   .mu_ground  Coulomb friction coefficient      [-]

if nargin < 5, recover = 0; end

% =========================================================================
% 1. PRE-PROCESSING
% =========================================================================
points_raw    = [rod_data(:,1:3); rod_data(:,4:6)];
tol_geo       = 1e-6;
points_rounded = round(points_raw / tol_geo) * tol_geo;
[unique_nodes, ~, ic] = unique(points_rounded, 'stable', 'rows');
num_nodes      = size(unique_nodes, 1);
num_beams_total = size(rod_data, 1);
is_rod         = rod_data(:,7) < 1;
num_actual_beams = sum(~is_rod);
node_ids_start = ic(1:num_beams_total);
node_ids_end   = ic(num_beams_total+1:end);

% Lumped masses
node_masses = zeros(num_nodes, 1);
for k = 1:num_beams_total
    if is_rod(k), prop_id = round(rod_data(k,7)*10); else, prop_id = rod_data(k,7); end
    p = properties(prop_id);
    L = norm(rod_data(k,4:6) - rod_data(k,1:3));
    m_half = (p.rho * p.A * L) / 2;
    node_masses(node_ids_start(k)) = node_masses(node_ids_start(k)) + m_half;
    node_masses(node_ids_end(k))   = node_masses(node_ids_end(k))   + m_half;
end

[~, top_idx]     = max(unique_nodes(:,3));
[~, z_sorted_idx]= sort(unique_nodes(:,3));
feet_indices     = z_sorted_idx(1:3);   % 3 lowest nodes = feet
p_cm             = unique_nodes(top_idx,:);
id_node_struct   = 1000 + top_idx;

S = Simparam;

% Absorber design
S_max      = 0.10;
F_crush_val = (0.5 * S.M_payload * S.mass_multiplier * S.v0(3)^2) / S_max;
K_stop_val  = 1e3;
V_lim_val   = 1e-1;
b_abs       = 567;   % Ns/m pure-viscous absorber

% Foot node MBDyn labels
foot_ids = 1000 + feet_indices;   % [f1 f2 f3]

% =========================================================================
% 2. WRITE FILE
% =========================================================================
fid = fopen(filename, 'w');

% --- header ---
fprintf(fid, 'begin: data;\n    problem: initial value;\nend: data;\n\n');

fprintf(fid, 'begin: initial value;\n');
fprintf(fid, '    initial time: %f;\n', S.t0);
fprintf(fid, '    final time:   %f;\n', S.tf);
fprintf(fid, '    time step:    %e;\n', S.dt);
fprintf(fid, '    method: ms, 0.6;\n');
fprintf(fid, '    linear solver: naive;\n');
fprintf(fid, '    max iterations: %d;\n', S.max_iterations);
fprintf(fid, '    tolerance: %e;\n', S.simtol);
fprintf(fid, '    derivatives coefficient: %e;\n', S.derivatives_coeff);
fprintf(fid, 'end: initial value;\n\n');

% --- scalar parameters ---
fprintf(fid, 'set: real mass_multiplier = %f;\n', S.mass_multiplier);
fprintf(fid, 'set: real x0 = %f; set: real y0 = %f; set: real z0 = %f;\n', S.p0(1), S.p0(2), S.p0(3));
fprintf(fid, 'set: real vx0 = %f; set: real vy0 = %f; set: real vz0 = %f;\n', S.v0(1), S.v0(2), S.v0(3));
fprintf(fid, 'set: real wx = %f; set: real wy = %f; set: real wz = %f;\n', S.w0(1), S.w0(2), S.w0(3));
fprintf(fid, 'set: real KDF = %f; set: real csi = %f;\n', S.KDF, S.csi);
fprintf(fid, 'set: real g = %f; set: const real LROD = 1.0;\n', S.g);

% payload inertia
fprintf(fid, 'set: real M_payload = %f * mass_multiplier;\n', S.M_payload);
fprintf(fid, 'set: const real R_pay = %f; set: const real H_pay = %f;\n', S.R_payload, S.H_payload);
fprintf(fid, 'set: real Jxx_pay = M_payload * (3.0*R_pay^2 + H_pay^2) / 12.0;\n');
fprintf(fid, 'set: real Jzz_pay = M_payload * R_pay^2 / 2.0;\n');

% absorber
fprintf(fid, 'set: const real S_max = %f; set: const real F_crush = %f;\n', S_max, F_crush_val);
fprintf(fid, 'set: const real K_stop = %e; set: const real V_lim = %f;\n', K_stop_val, V_lim_val);
fprintf(fid, 'set: const real b = %e;\n\n', b_abs);

% Bekker soil parameters (printed as literals into the element card)
kc        = S.kc;
kphi      = S.kphi;
n_bek     = S.n_bek;
b_pad     = S.b_pad;
k_plastic = S.k_plastic;
h_hydro   = S.h_hydro;
VN0       = S.VN0;
ZN0       = S.ZN0;
mu_s      = S.mu_ground;

% --- control data ---
% nodes  : 1 (ground 99) + num_nodes (struct 1000+) + 1 (payload 5000)
% bodies : num_nodes (struct) + 1 (payload)
% joints : sum(is_rod) rod joints + 1 clamp(99) + 1 prismatic absorber
% forces : 1 (absorber follower)
% loadable: 1 (surface impact)
fprintf(fid, 'begin: control data;\n');
fprintf(fid, '    output frequency: %d;\n', S.out_freq);
fprintf(fid, '    structural nodes: %d;\n', num_nodes + 2);
fprintf(fid, '    rigid bodies: %d;\n',     num_nodes + 1);
fprintf(fid, '    joints: %d;\n',           sum(is_rod) + 2);
fprintf(fid, '    beams: %d;\n',            num_actual_beams);
fprintf(fid, '    forces: 1;\n');
fprintf(fid, '    loadable elements: 1;\n');
fprintf(fid, '    gravity;\n');
fprintf(fid, '    module load: "libmodule-soil-plasticity";\n');
fprintf(fid, '    default output: structural, accelerations, loadable;\n');
fprintf(fid, 'end: control data;\n\n');

% --- nodes ---
fprintf(fid, 'begin: nodes;\n');
fprintf(fid, '    structural: 99, static, 0., 0., -LROD, eye, null, null;\n');
for i = 1:num_nodes
    dx = unique_nodes(i,1)-p_cm(1);
    dy = unique_nodes(i,2)-p_cm(2);
    dz = unique_nodes(i,3)-p_cm(3);
    fprintf(fid, '    structural: %d, dynamic, %f+x0, %f+y0, %f+z0, eye,\n', ...
        1000+i, unique_nodes(i,1), unique_nodes(i,2), unique_nodes(i,3));
    fprintf(fid, '        vx0+wy*%f-wz*%f, vy0+wz*%f-wx*%f, vz0+wx*%f-wy*%f, wx, wy, wz;\n', ...
        dz, dy, dx, dz, dy, dx);
end
% payload node
fprintf(fid, '    structural: 5000, dynamic, %f+x0, %f+y0, %f+z0+H_pay/2., eye,\n', p_cm(1), p_cm(2), p_cm(3));
fprintf(fid, '        vx0, vy0, vz0, wx, wy, wz;\n');
fprintf(fid, 'end: nodes;\n\n');

% --- plugin variables (needed only for absorber) ---
fprintf(fid, 'set: [node, z_struct, %d, structural, string="X[3]"];\n',  id_node_struct);
fprintf(fid, 'set: [node, vz_struct, %d, structural, string="XP[3]"];\n', id_node_struct);
fprintf(fid, 'set: [node, z_pay,  5000, structural, string="X[3]"];\n');
fprintf(fid, 'set: [node, vz_pay, 5000, structural, string="XP[3]"];\n');
fprintf(fid, '\n');

% --- elements ---
fprintf(fid, 'begin: elements;\n');
fprintf(fid, '    gravity: uniform, 0., 0., -1., const, g;\n');
fprintf(fid, '    joint: 99, clamp, 99, node, eye;\n\n');

% structural bodies
for i = 1:num_nodes
    m_b = node_masses(i) * S.mass_multiplier;
    fprintf(fid, '    body: %d, %d, %f, null, diag, 0.1, 0.1, 0.1;\n', 3000+i, 1000+i, m_b);
end

% payload body
fprintf(fid, '    body: 3999, 5000, M_payload, 0., 0., 0.,\n');
fprintf(fid, '        diag, Jxx_pay, Jxx_pay, Jzz_pay;\n\n');

% absorber between top structural node and payload
fprintf(fid, '    joint: 6000, prismatic, %d, 5000;\n', id_node_struct);
fprintf(fid, '    force: 6001, follower internal, %d, position, null, 5000, position, null, 0.,0.,-1.,\n', id_node_struct);
switch S.absorber_type
    case 0  % plastic crusher
        fprintf(fid, '        string, "max(0.,tanh((vz_struct-vz_pay)/V_lim))*(F_crush+(z_struct-z_pay>S_max)*K_stop*(z_struct-z_pay-S_max))+(z_struct-z_pay<0.)*K_stop*(z_struct-z_pay)";\n');
    case 1  % pure viscous
        fprintf(fid, '        string, "max(0.,tanh((vz_struct-vz_pay)/V_lim))*(vz_struct-vz_pay)*b";\n');
    case 2  % spring-damper
        fprintf(fid, '        string, "max(0.,tanh((vz_struct-vz_pay)/V_lim))*((vz_struct-vz_pay)*b+K_stop*(z_struct-z_pay))";\n');
    case 3  % viscous + friction pad
        fprintf(fid, '        string, "max(0.,tanh((vz_struct-vz_pay)/V_lim))*(vz_struct-vz_pay)*b+40.*min(1.,(vz_struct-vz_pay)/%e)";\n', S.vtol);
end
fprintf(fid, '\n');

% --- surface impact element (replaces all shadow/extra node contact machinery) ---
fprintf(fid, '    # Bekker-Wong soil contact on 3 feet\n');
fprintf(fid, '    # Params: kc=%.1f  kphi=%.1f  n=%.2f  b=%.4f m\n', kc, kphi, n_bek, b_pad);
fprintf(fid, '    # k_plastic=%.2f  h=%.1f  VN0=%.3f  ZN0=%.4f  mu=%.3f\n', k_plastic, h_hydro, VN0, ZN0, mu_s);
fprintf(fid, '    user defined: 9001, surface impact,\n');
fprintf(fid, '        plane normal, 0., 0., 1.,\n');
fprintf(fid, '        plane point,  0., 0., 0.,\n');
fprintf(fid, '        kc,           %.6g,\n', kc);
fprintf(fid, '        kphi,         %.6g,\n', kphi);
fprintf(fid, '        n,            %.6g,\n', n_bek);
fprintf(fid, '        pad radius,   %.6g,\n', b_pad);
fprintf(fid, '        plastic ratio, %.6g,\n', k_plastic);
if h_hydro > 0
    fprintf(fid, '        hydrodynamic,\n');
    fprintf(fid, '            h,   %.6g,\n', h_hydro);
    fprintf(fid, '            VN0, %.6g,\n', VN0);
    fprintf(fid, '            ZN0, %.6g,\n', ZN0);
end
fprintf(fid, '        friction, %.6g,\n', mu_s);
fprintf(fid, '        nodes, 3,\n');
fprintf(fid, '            %d, %d, %d;\n\n', foot_ids(1), foot_ids(2), foot_ids(3));

% --- structural elements (rods and beams) ---
for k = 1:num_beams_total
    n1 = 1000 + node_ids_start(k);
    n2 = 1000 + node_ids_end(k);
    P1 = unique_nodes(node_ids_start(k),:);
    P2 = unique_nodes(node_ids_end(k),:);
    L_elem = norm(P2 - P1);
    L_ref  = rod_data(k, 8);

    if is_rod(k)
        p      = properties(round(rod_data(k,7)*10));
        EA     = p.E * p.A;
        rhoA   = p.rho * p.A;
        omega_axial = pi * sqrt(EA / (rhoA * L_ref^2));
        eta_axial   = 2 / omega_axial;
        c_axial_sect = eta_axial * EA;
        fprintf(fid, '    joint: %d, rod, %d, %d, from nodes,\n', 2000+k, n1, n2);
        fprintf(fid, '        linear viscoelastic, KDF*%e, csi*%e;\n', EA, c_axial_sect);
    else
        p     = properties(rod_data(k,7));
        delta = P2 - P1;
        e1    = delta' / L_elem;
        perp  = null(e1');
        e2    = perp(:,1); e3 = perp(:,2);
        R     = [e1, e2, e3];

        EA = p.E*p.A; GA = p.G*p.A; GJ = p.G*p.J; EI = p.E*p.I;
        rhoA = p.rho*p.A; rhoI = p.rho*p.I; rhoJ = p.rho*p.J;

        omega_axial   = pi*sqrt(EA/(rhoA*L_ref^2));
        omega_shear   = pi*sqrt(GA/(rhoA*L_ref^2));
        omega_torsion = pi*sqrt(GJ/(rhoJ*L_ref^2));
        beta1         = 0.597*pi;
        omega_bending = beta1^2*sqrt(EI/(rhoA*L_ref^4));

        c_ax  = (2*S.csi/omega_axial)   * EA;
        c_sh  = (2*S.csi/omega_shear)   * GA;
        c_tor = (2*S.csi/omega_torsion)  * GJ;
        c_ben = (2*S.csi/omega_bending)  * EI;

        fprintf(fid, '    beam2: %d, %d, null, %d, null, matr,\n', 2000+k, n1, n2);
        fprintf(fid, '        %.16e, %.16e, %.16e,\n', R(1,1), R(1,2), R(1,3));
        fprintf(fid, '        %.16e, %.16e, %.16e,\n', R(2,1), R(2,2), R(2,3));
        fprintf(fid, '        %.16e, %.16e, %.16e,\n', R(3,1), R(3,2), R(3,3));
        fprintf(fid, '        linear viscoelastic generic,\n');
        fprintf(fid, '        diag, KDF*%e, KDF*%e, KDF*%e, KDF*%e, KDF*%e, KDF*%e,\n', EA,GA,GA,GJ,EI,EI);
        fprintf(fid, '        diag, csi*%e, csi*%e, csi*%e, csi*%e, csi*%e, csi*%e;\n', c_ax,c_sh,c_sh,c_tor,c_ben,c_ben);
    end
end

fprintf(fid, 'end: elements;\n\n');
fclose(fid);
end
