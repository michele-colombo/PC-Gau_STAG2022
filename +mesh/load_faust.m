function M = load_faust(mesh_name, n_eigs, shuffle)

    if nargin < 3
        shuffle = false;
    end
    if nargin < 2
        n_eigs = 200;
    end

    mesh_dir = fullfile('Meshes', 'FAUST');
    M = mesh.init(fullfile(mesh_dir, [mesh_name,'.mat']));
    if shuffle        
        M = mesh.transform.shuffle_shape(M);
    else
        M.gt = [1:M.n]';
        M.inv_gt = [1:M.n]';
    end
        
    %Src = mesh.transform.rotate(Src, 'x', 90);      
    %Src = mesh.transform.rotate(Src, 'y', -40);
    M = mesh.transform.normalize(M);    
    [M.S, ~, M.A] = mesh.proc.FEM_higher(M, 1, 'Dirichlet');
    [M.Phi, M.Lambda] = eigs(M.S, M.A, n_eigs, -1e-5);
    M.Lambda = diag(M.Lambda);
    M.evals = M.Lambda;
    M.evecs = M.Phi;
    
end
