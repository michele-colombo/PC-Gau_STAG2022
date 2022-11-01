clear all; close all; clc; rng(0); init_path;

%% parameters
M_name = 'tr_reg_023';

n_points = 1000;    % number of vertices used to build PC-Gau (q)
sigma = 0.05;       % amplitude of the gaussians used to build PC-Gau (sigma)
k = 60;             % number of atoms to consider

s = 80;            % sets the scale for the egdc metric
t = 10;             % sets the scale for the mgd metric

%% load meshes and compute LB basis
M = mesh.load_faust(M_name);

%% compute PC-Gau
para.n_atoms = 200;
M.our_basis = compute_PC_Gau(M, n_points, sigma, true, para);

%% embeddings
% embeddings with delta functions are the basis matrices transposed
embeddings_ours = M.our_basis(:,1:k)';
embeddings_LB = M.Phi(:,1:k)';

% For each vertex, compute the nearest embedding and its embedding distance
[min_dist_ours, ours_nn] = compute_min_dist(embeddings_ours);    
[min_dist_LB, LB_nn] = compute_min_dist(embeddings_LB);   

% initialize vectors to contain metric values
Dis_ours = zeros(M.n, 1);
Dis_LB = zeros(M.n, 1);
EGDC_ours = zeros(M.n, 1);
EGDC_LB = zeros(M.n, 1);
MGD_ours = zeros(M.n, 1);
MGD_LB = zeros(M.n, 1);
        
%% Discrimination power (Dis)
disp('computing metrics, it may take a while...');
for point = 1:M.n
    % geodesic distances from the current point
    geodesic_distances = distances(M.norm_G, point)';  
    % embedding distances from the current point, for ours and LB
    emb_dist_ours = pdist2(embeddings_ours(:,point)', embeddings_ours')';
    emb_dist_LB = pdist2(embeddings_LB(:,point)', embeddings_LB')';
    
    % Set the metrics' values in the current point
    
    % Discrimination power: distance from the nearest embedding normalized
    % by the geodesic distance of the corresponding vertex
    Dis_ours(point) = min_dist_ours(point) / geodesic_distances(ours_nn(point));
    Dis_LB(point) = min_dist_LB(point) / geodesic_distances(LB_nn(point));
    
    % EGDC: correlation between geodesic and embedding distances, for the s
    % nearest embeddings
    EGDC_ours(point) = compute_egdc(geodesic_distances, emb_dist_ours, s);
    EGDC_LB(point) = compute_egdc(geodesic_distances, emb_dist_LB, s);
    
    % MGD: mean geodesic distance of the t nearest embeddings, normalized
    % by the mean geodesic distance of the actual t nearest vertices
    MGD_ours(point) = compute_mgd(geodesic_distances, emb_dist_ours, t);
    MGD_LB(point) = compute_mgd(geodesic_distances, emb_dist_LB, t);   
end

%% figure parameters
viz_params.x = 83;
viz_params.y = 0;
viz_params.z = -10;
viz_params.view = [-37.5 10];
viz_params.diffuseStrength = 0.5;

%% Plot Dis
fig_Dis = figure('position', [100 100 900 900]);
lim_inf = min([Dis_ours; Dis_LB]);
lim_sup = max([Dis_ours; Dis_LB])*0.5;      % rescaled for better visualization

subplot(1,2,1); utils.plot_scalar_map_2(M, Dis_ours, hot, viz_params);
caxis([lim_inf, lim_sup]); 
colorbar('location', 'southoutside');
title(['Dis (ours), mean: ', num2str(mean(Dis_ours), '%.1f')]);

subplot(1,2,2); utils.plot_scalar_map_2(M, Dis_LB, hot, viz_params);
caxis([lim_inf, lim_sup]); 
colorbar('location', 'southoutside');
title(['Dis (LB), mean: ', num2str(mean(Dis_LB), '%.1f')]);

sgtitle({'\fontsize{22}Discrimination power'; '\fontsize{16}Lower (darker) is worse'})

%% Plot EGDC
fig_EGDC = figure('position', [100 100 900 900]);
lim_inf = min([EGDC_ours; EGDC_LB]);
lim_sup = max([EGDC_ours; EGDC_LB]);

subplot(1,2,1); utils.plot_scalar_map_2(M, EGDC_ours, hot, viz_params);
caxis([lim_inf, lim_sup]);
colorbar('location', 'southoutside'); 
title(['EGDC (ours), mean: ', num2str(mean(EGDC_ours), '%.3f')]);

subplot(1,2,2); utils.plot_scalar_map_2(M, EGDC_LB, hot, viz_params);
caxis([lim_inf, lim_sup]); 
colorbar('location', 'southoutside');
title(['EGDC (LB), mean: ', num2str(mean(EGDC_LB), '%.3f')]);

sgtitle({'\fontsize{22}Embedding/Geodeisc distance correlation'; '\fontsize{16}Lower (darker) is worse'})

%% Plot MGD
fig_MGD = figure('position', [100 100 900 900]);
lim_inf = min([MGD_ours; MGD_LB]);
lim_sup = max([MGD_ours; MGD_LB])*0.5;      %rescaled for better visualization

subplot(1,2,1); utils.plot_scalar_map_2(M, MGD_ours, flipud(hot), viz_params);
caxis([lim_inf, lim_sup]); 
colorbar('location', 'southoutside');
title(['MGD (ours), mean: ', num2str(mean(MGD_ours), '%.2f')]);

subplot(1,2,2); utils.plot_scalar_map_2(M, MGD_LB, flipud(hot), viz_params);
caxis([lim_inf, lim_sup]); 
colorbar('location', 'southoutside');
title(['MGD (LB), mean: ', num2str(mean(MGD_LB), '%.2f')]);

sgtitle({'\fontsize{22}Mean Geodesic Distance (normalized)'; '\fontsize{16}Higher (darker) is worse'})



%%
function egdc = compute_egdc(geodesic_distances, spectral_distances, s)
    [sorted_spectral, idxs] = sort(spectral_distances);
    sorted_geodesic = geodesic_distances(idxs);
    egdc = corr(sorted_spectral(1:s), sorted_geodesic(1:s));

end

function mgd_norm = compute_mgd(geodesic_distances, spectral_distances, t)
    [~, idxs] = sort(spectral_distances);
    sorted_geodesic = geodesic_distances(idxs);
    real_geo_dist = mink(geodesic_distances, t);
    mgd = mean(sorted_geodesic(1:t));
    mgd_norm = mgd ./ mean(real_geo_dist);
end

function [min_dist, nn] = compute_min_dist(embeddings)
    [nn,min_dist] = knnsearch(embeddings', embeddings', 'K', 2);
    min_dist = min_dist(:,2);
    nn = nn(:,2);
end
