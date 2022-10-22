function M = load_mat(filename)

    tmp = load(filename);
    if isfield(tmp, 'surface')
        shape = tmp.surface;
    elseif isfield(tmp, 'shape')
        shape = tmp.shape;
    elseif isfield(tmp, 'Shape_df')
        shape = tmp.Shape_df;
    end
    

    
    if isfield(shape, 'VERT')
        M.VERT = shape.VERT;
    elseif isfield(shape, 'X')
        M.VERT = [shape.X, shape.Y, shape.Z];
    end
    
    if isfield(shape, 'corr')
        M.corr = shape.corr;
    elseif isfield(shape, 'corr')
        M.corr = tmp.corr;
    end
        
    M.TRIV = shape.TRIV; 
    M.n = size(M.VERT, 1);
    M.m = size(M.TRIV, 1);
    
    if isfield(tmp, 'G')
        M.G = tmp.G;
    end
    
    if isfield(tmp, 'norm_G')
        M.norm_G = tmp.norm_G;
    end    
    
    if isfield(tmp, 'name')
        M.name = tmp.name;
    end
    
    if isfield(tmp, 'smpl_matches')
        M.smpl_matches = tmp.smpl_matches;
    end  
    
end