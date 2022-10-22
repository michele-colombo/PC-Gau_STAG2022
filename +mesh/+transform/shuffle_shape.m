function S = shuffle_shape(M)

n = size(M.VERT,1);
S = {};

perm = randperm(n);
gt(perm) = 1:n;
S.VERT = M.VERT(perm,:);

if isfield(M, 'TRIV')
    S.TRIV = gt(M.TRIV);
    S.m = size(S.TRIV,1);
elseif isfield(M, 'EDGES')
    S.EDGES = gt(M.EDGES);
    S.m = size(S.EDGES,1);
else
    error('Shape must have either triangular or edge faces.')
end

[~,S.gt] = ismember(perm,1:n);

S.n = size(S.VERT,1);
S.inv_gt = gt;

end
