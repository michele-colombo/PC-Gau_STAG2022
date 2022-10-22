function M = load_shrec20(mesh_name, n_eigs)
    if nargin < 2
        n_eigs = 200;
    end

    M = mesh.init(fullfile('Meshes', 'SHREC20b_lores', 'normalized_meshes', [mesh_name,'.mat']));
    M = mesh.transform.normalize(M);    
    [M.S, ~, M.A] = mesh.proc.FEM_higher(M, 1, 'Neumann');
    [M.Phi, M.Lambda] = eigs(M.S, M.A, n_eigs, -1e-5);
    M.Lambda = diag(M.Lambda);
    M.evals = M.Lambda;
    M.evecs = M.Phi;
    
    gt_file = fullfile('Meshes','SHREC20b_lores', 'corres', [mesh_name,'.mat']);
    M.GT = load(gt_file);
    
%     sf = mesh.proc.compute_indicatorBasis(M, 1e-5, 1, M.GT.verts, [], false, n_eigs);
%     display_basis(M, sf, num2cell(M.GT.points), size(sf, 2));
end