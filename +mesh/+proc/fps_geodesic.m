function [ S ] = fps_geodesic(M, K, seed)
% FPS_GEODESIC(M, K, seed) samples K vertices from the vertices of M by using farthest point sampling.
% The farthest point sampling starts with vertex seed and uses the geodesic
% distance between the vertices.
% -  M is a mesh with a field named 'G' or 'norm_G' representing its mesh
%    graph
% -  K is the number of samples
% -  seed is the index of the first vertex in the sample set. (1<=seed<=K)
% Returns
% -  S is a K-dimensional vector that includes the indeces of the K sample
%    vertices.

S = zeros(K,1);
S(1) = seed;
if isfield(M, 'norm_G')
    mesh_graph = M.norm_G;
elseif isfield(M, 'G')
    mesh_graph = M.G;
else
    disp('The provided mesh has no mesh graph field');
    return
end

d = distances(mesh_graph, seed)';

for i=2:K
    [~,m] = max(d);
    S(i) = m(1);
    d = min(distances(mesh_graph, S(i))',d);
end

end