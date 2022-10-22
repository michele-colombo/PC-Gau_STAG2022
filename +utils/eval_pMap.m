function err = eval_pMap(S1, S2, T12, type, gt_corres)
% compute the error of the pMap
% Input:
%   T12: S1 -> S2
%   gt_corres: m-by-2 matrix, the m corresponding vertices (ground-truth to measure the T12)
%   type: 'geodesic' or ';euclidean' distance
if nargin < 4, type = 'euclidean'; end;

if nargin < 5
    gt_corres = [(1:S1.n)', (1:S1.n)'];
end

T12 = reshape(T12,[],1);
vid1 = gt_corres(:,1); % vtx on S1
vid2 = gt_corres(:,2); % vtx on S2 corresponding to vid1

if(strcmp(type, 'geodesic'))
    if isfield(S2,'Gamma') % we have geodesic distance matrix
        get_mat_entry = @(M,I,J) M(sub2ind(size(M),I,J));
        dists = get_mat_entry(S2.Gamma,vid2,T12(vid1));
    else
        dists = geodesics_pairs(S2, [T12(vid1), vid2]);
    end
elseif(strcmp(type, 'euclidean'))
    D = [S2.X S2.Y S2.Z];
    dists = sqrt(sum((D(vid2, :) - D(T12(vid1), :)).^2,2));
elseif(strcmpi(type, 'dijkstra'))
    if isfield(S2,'Gamma') % we have dijkstra distance matrix
        get_mat_entry = @(M,I,J) M(sub2ind(size(M),I,J));
        dists = get_mat_entry(S2.Gamma,vid2,T12(vid1));
    elseif isfield(S2, 'norm_G')
        c_matches = T12(vid1);
        c_gt_matches = vid2;
        [src_idx, src_ia, src_ic]=unique(c_matches);
        [sink_idx, sink_ia, sink_ic]=unique(c_gt_matches);
        geo_d = distances(S2.norm_G,src_idx,sink_idx);
        dist_matrix = geo_d(src_ic,sink_ic);
        [ii, jj] = find(isinf(dist_matrix));
        if not(isempty(ii))
            eucl=vecnorm(S2.VERT(c_matches(ii),:)-S2.VERT(c_gt_matches(jj),:),2,2);
            idx = sub2ind(size(dist_matrix),ii,jj);
            dist_matrix(idx)=eucl;
        end
        dists = diag(dist_matrix);
    else
        dists = dijkstra_pairs(S2, [T12(vid1), vid2]);
    end
else
    error('unknown type %s\n',type);
end

err = dists/sqrt(sum(mesh.proc.tri_areas(S2)));

end
