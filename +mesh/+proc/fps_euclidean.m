function [ S ] = fps_euclidean(V, K, seed)
% FPS_EUCLIDEAN(V, K, seed) samples K vertices from V by using farthest point sampling.
% The farthest point sampling starts with vertex seed and uses the euclidean
% metric of the 3-dimensional embedding space.
% -  V is a n-by-3 matrix storing the positions of n vertices
% -  K is the number of samples
% -  seed is the index of the first vertex in the sample set. (1<=seed<=K)
% Returns
% -  S is a K-dimensional vector that includes the indeces of the K sample
%    vertices.

S = zeros(K,1);
S(1) = seed;
d = pdist2(V,V(seed,:));

for i=2:K
    [~,m] = max(d);
    S(i) = m(1);
    d = min(pdist2(V,V(S(i),:)),d);
end

end