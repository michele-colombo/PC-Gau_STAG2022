clear all; close all; clc; rng(0); init_path;

%% parameters
M_name = 'tr_reg_007';
N_name = 'tr_reg_020';
C_computation_method = 'GT';      % 'GT' to compute the functional map C from groundtruth, 
                                  % 'NO17' to compute C using NO17

n_points = 1000;
sigma = 0.05;
k = 60;


%% load meshes and compute LB basis
M = mesh.load_faust(M_name);
N = mesh.load_faust(N_name);


%% compute PC-Gau
para.n_atoms = 200;
M.our_basis = compute_PC_Gau(M, n_points, sigma, true, para);
N.our_basis = compute_PC_Gau(N, n_points, sigma, true, para);


%% visualize basis
n_tiles = 8;
fig_basis = figure('position', [100 100 1800 900]);
tlo = tiledlayout(1,n_tiles);
for i=1:n_tiles
    nexttile(tlo);
    utils.plot_scalar_map(M, M.our_basis(:,i), utils.cmaps.bwr);
    title(['atom ', num2str(i)]);
end
title(tlo, 'PC-Gau basis', ' ', 'fontSize', 30)


%% compute C (with ground-truth correspondence)
T_corres = [N.inv_gt, M.inv_gt];                     
LM = (500:1000:6000)';

if strcmpi(C_computation_method, 'NO17')
    % compute functional map C with NO17, for ours and LB
    C_ours = build_funmap_basis(M,N,M.our_basis(:,1:k), N.our_basis(:,1:k), 200, 200, LM, LM);
    C_LB = build_funmap_basis(M,N,M.Phi(:,1:k), N.Phi(:,1:k), 200, 200, LM, LM);
else
    % compute functional map C using the ground truth, for ours and LB
    C_ours = N.our_basis(T_corres(:,1), 1:k)\M.our_basis(T_corres(:,2), 1:k);       
    C_LB = N.Phi(T_corres(:,1), 1:k)\M.Phi(T_corres(:,2), 1:k);                            
end


%% convert fMap C to pMap T
T_ours = utils.fMap2pMap(M.our_basis(:,1:k), N.our_basis(:,1:k), C_ours);
T_LB = utils.fMap2pMap(M.Phi(:,1:k), N.Phi(:,1:k), C_LB);

error_ours = utils.eval_pMap(N,M,T_ours, 'dijkstra', T_corres);
age_ours = mean(error_ours);
error_LB = utils.eval_pMap(N,M,T_LB, 'dijkstra', T_corres);
age_LB = mean(error_LB);


%% plot map and errors
lim_inf = min([error_ours; error_LB]);
lim_sup = max([error_ours; error_LB]);
viz_params.x = 83;
viz_params.y = 0;
viz_params.z = -10;
viz_params.view = [-37.5 10];
viz_params.diffuseStrength = 0.5;

fig_matching = figure('position', [100 100 1800 700]);
tlo = tiledlayout(1,5);

nexttile(); 
utils.visualize_map_on_source(N, M, viz_params);
title('source');

nexttile(); utils.visualize_map_on_target(N, M, T_ours, viz_params);
title('pMap (ours)');

nexttile(); utils.visualize_map_on_target(N, M, T_LB, viz_params);
title('pMap (LB)');

nexttile(); utils.plot_scalar_map_2(N, error_ours, flipud(hot), viz_params);
caxis([lim_inf, lim_sup]); 
title(['error (ours): ', num2str(age_ours*1000, '%.1f'), '\cdot 10^{-3}']);

nexttile(); utils.plot_scalar_map_2(N, error_LB, flipud(hot), viz_params);
caxis([lim_inf, lim_sup]); 
colorbar;
title(['error (LB): ', num2str(age_LB*1000, '%.1f'), '\cdot 10^{-3}']);

title(tlo, 'Example of shape matching using PC-Gau and LB', 'fontSize', 24);
subtitle(tlo, ['C computed with ', C_computation_method], 'fontSize', 18);



