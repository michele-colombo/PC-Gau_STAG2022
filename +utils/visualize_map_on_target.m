function colors = visualize_map_on_target(S1, S2, T12, para)
if (nargin >= 4 && isfield(para, 'x'))
    S1 = mesh.transform.rotate(S1, 'x', para.x);    
    S2 = mesh.transform.rotate(S2, 'x', para.x);
end
if (nargin >= 4 && isfield(para, 'y'))
    S1 = mesh.transform.rotate(S1, 'y', para.y);    
    S2 = mesh.transform.rotate(S2, 'y', para.y);
end
if (nargin >= 4 && isfield(para, 'z'))
    S1 = mesh.transform.rotate(S1, 'z', para.z);    
    S2 = mesh.transform.rotate(S2, 'z', para.z);
end

[g1,g2,g3] = set_mesh_color(S2);
    f1 = g1(T12);
    f2 = g2(T12);
    f3 = g3(T12);
    colors = [f1, f2, f3];
    trimesh(S1.TRIV, S1.X, S1.Y, S1.Z, ...
    'FaceVertexCData', colors,...
    'FaceColor', 'interp', 'EdgeColor', 'none'); axis equal; axis off; 
    if nargin>=4 && isfield(para, 'view')
        view(para.view);
    end    
end


function [g1,g2,g3] = set_mesh_color(S)
g1 = normalize_function(0,1,S.X);
g2 = normalize_function(0,1,S.Y);
g3 = normalize_function(0,1,S.Z);
g1 = reshape(g1,[],1);
g2 = reshape(g2,[],1);
g3 = reshape(g3,[],1);
end

function fnew = normalize_function(min_new,max_new,f)
    fnew = f - min(f);
    fnew = (max_new-min_new)*fnew/max(fnew) + min_new;
end