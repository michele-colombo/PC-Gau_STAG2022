function M = load_faust_remeshed(mesh_name, n_eigs)
    if nargin < 2
        n_eigs = 200;
    end

    mesh_dir = fullfile('Meshes', 'FAUST_remeshed');
    M = mesh.init(fullfile(mesh_dir, [mesh_name,'.mat']));
    M.inv_gt = importdata(fullfile(mesh_dir, [mesh_name,'.vts']));
    [~,M.gt] = ismember(1:M.n, M.inv_gt);
    M.gt = M.gt';
    

    M = mesh.transform.normalize(M);    
    [M.S, ~, M.A] = mesh.proc.FEM_higher(M, 1, 'Dirichlet');
    [M.Phi, M.Lambda] = eigs(M.S, M.A, n_eigs, -1e-5);
    M.Lambda = diag(M.Lambda);
    M.evals = M.Lambda;
    M.evecs = M.Phi;
    
end