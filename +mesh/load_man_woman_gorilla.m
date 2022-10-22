function M = load_man_woman_gorilla(mesh_name, n_eigs)
    if nargin < 2
        n_eigs = 200;
    end

    M = mesh.init(fullfile('Meshes', 'man_woman_gorilla', 'normalized_meshes', [mesh_name,'.mat']));
%     load(fullfile('Meshes', 'man_woman_gorilla', [mesh_name,'.mat']));
    M.inv_gt = M.corr;
    [~,M.gt] = ismember(1:M.n, M.inv_gt);
    M.gt = M.gt';

    M = mesh.transform.normalize(M);    
    [M.S, ~, M.A] = mesh.proc.FEM_higher(M, 1, 'Neumann');
    [M.Phi, M.Lambda] = eigs(M.S, M.A, n_eigs, -1e-5);
    M.Lambda = diag(M.Lambda);
    M.evals = M.Lambda;
    M.evecs = M.Phi;

end