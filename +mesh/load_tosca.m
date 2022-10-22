function M = load_tosca(mesh_name, n_eigs, shuffle, resolution)
    if nargin < 2
        n_eigs = 200;
    end
    if nargin < 4
        resolution = 'high';
    end

    if strcmpi(resolution, 'low')
        M = mesh.init(fullfile('Meshes', 'TOSCA_lores', [mesh_name,'.mat']));
    elseif strcmpi (resolution, 'high')
        M = mesh.init(fullfile('Meshes', 'TOSCA_hires', 'normalized_meshes', [mesh_name,'.mat']));
    end
    
    %load(fullfile('Meshes', 'man_woman_gorilla', [mesh_name,'.mat']));
    M.gt = [1:M.n]';
    M.inv_gt = [1:M.n]';

    M = mesh.transform.normalize(M);    
    [M.S, ~, M.A] = mesh.proc.FEM_higher(M, 1, 'Neumann');
    [M.Phi, M.Lambda] = eigs(M.S, M.A, n_eigs, -1e-5);
    M.Lambda = diag(M.Lambda);
    M.evals = M.Lambda;
    M.evecs = M.Phi;

end