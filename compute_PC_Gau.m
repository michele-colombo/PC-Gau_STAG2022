function [P, selected_points] = compute_PC_Gau(M, points, sigma, normalize, para, D)
%COMPUTE_PC_GAU compute the PC-Gau basis for a given mesh M
% - M is the input mesh with M.n vertices
% - points is either the number of vertices to sample for the basis
% computation or a vector of the indices of the vertices to use
% - sigma determines the amplitude of the gaussians. It can be either a
% scalar (it uses the same sigma for all vertices) or a size M.n vector of
% sigma values (one for each vertex)
% - normalize is a flag that enables/disables normalization of the gaussian
% prior to PCA computation
% - para contains other parameters of the computation, such as:
%    - point_selection_method: determines the methods used to sample the
%    vertices. It can be 'euclidean', 'random' or 'geodesic'. This is
%    meaningful only when points is a scalar
%    - n_atoms: determines the number of atoms produced in output
% - P is a M.n-by-k matrix containing the basis atoms as column vectors. If
% para.n_atoms has been specified, k<=n_atoms.
% - selected_points contains the indices of the vertices used to generate
% the gaussian functions
if nargin < 4
    normalize = true;
end

if nargin >= 5 && isfield(para, 'point_selection_method')
    point_selection_method = para.point_selection_method;
else
    point_selection_method = 'euclidean';
end


if isscalar(points)     % if points is a scalar, it is the number of points to select
    if strcmpi(point_selection_method, 'euclidean')
        selected_points = mesh.proc.fps_euclidean(M.VERT, points, 1);
    elseif strcmpi(point_selection_method, 'random')
        selected_points = randperm(M.n, points);
    elseif strcmpi(point_selection_method, 'geodesic')
        selected_points = mesh.proc.fps_geodesic(M, points, 1);
    end
else        % otherwise points contains the selected points
    selected_points = points;
end

if nargin < 6 || isempty(D)
    fprintf('computing distances...'); tic;
    if isfield(M, 'norm_G')
        D = distances(M.norm_G, selected_points)';
    elseif isfield(M, 'G')
        D = distances(M.G, selected_points)';
        D = D ./ sqrt(sum(mesh.proc.tri_areas(M)));
    else
        D = dijkstra_all_to_all(M, selected_points);
    end
    toc;
else
    D = D(:, selected_points);
end

if ~isscalar(sigma)     % if sigma is a vector, it represents the value of sigma to use in each vertex of the mesh
    if size(sigma,1)>size(sigma,2)
        sigma = sigma';
    end
    sigma = sigma(selected_points); 
end
gaussians = exp(-1 * (D.^2 ./ sigma));  % this handles both the case of a constant sigma (scalar) and of a vector of sigma values


if normalize
    norms = sqrt(diag(gaussians'*M.A*gaussians))';
    gaussians = gaussians ./ repmat(norms, [M.n 1]);
end


% fprintf('computing PCA... ');tic;
P = pca(gaussians', 'Centered', false, 'VariableWeights', full(diag(M.A)));
% toc;


if nargin >= 5 && isfield(para, 'n_atoms')     % outputs only the first n_atoms atoms
    if para.n_atoms<size(P,2)
        P = P(:,1:para.n_atoms);
    end
end


end