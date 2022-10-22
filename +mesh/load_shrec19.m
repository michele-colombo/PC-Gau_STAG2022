function M = load_shrec19(mesh_name, n_eigs, shuffle)

    if nargin < 3
        shuffle = false;
    end
    if nargin < 2
        n_eigs = 200;
    end

    M = mesh.init(fullfile('Meshes', 'SHREC19', 'normalized_meshes', [mesh_name,'.mat']));
    
    M.inv_gt = M.corr;
    if size(M.inv_gt,2)>size(M.inv_gt,1)
        M.inv_gt = M.inv_gt';
    end
    [~,M.gt] = ismember(1:M.n, M.inv_gt);
    M.gt = M.gt';

    M = mesh.transform.normalize(M);    
    [M.S, ~, M.A] = mesh.proc.FEM_higher(M, 1, 'Neumann');
    [M.Phi, M.Lambda] = eigs(M.S, M.A, n_eigs, -1e-5);
    M.Lambda = diag(M.Lambda);
    M.evals = M.Lambda;
    M.evecs = M.Phi;
    
end