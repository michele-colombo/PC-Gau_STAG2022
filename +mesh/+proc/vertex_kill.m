function N = vertex_kill(M, perc)

    N = M;
    
    while N.n > perc * M.n
        E = [N.TRIV(:, [1 2]); N.TRIV(:, [2 3]); N.TRIV(:, [1 3])];
        E = sort(E, 2);
        E = unique(E, 'rows');
        EVec = N.VERT(E(:, 1), :) - N.VERT(E(:, 2), :);
        ELen = vecnorm(EVec, 2);
        
        [~, i] = min(ELen);
        todel = E(i, :);
        NewVert = (N.VERT(todel(1), :) + N.VERT(todel(2), :)) / 2;
        N.TRIV(N.TRIV == todel(1)) = N.n + 1;
        N.TRIV(N.TRIV == todel(2)) = N.n + 1;
        N.TRIV(N.TRIV >= todel(2)) = N.TRIV(N.TRIV >= todel(2)) - 1;
        N.TRIV(N.TRIV >= todel(1)) = N.TRIV(N.TRIV >= todel(1)) - 1;
        N.VERT = [N.VERT(1:todel(1)-1, :);
                  N.VERT(todel(1)+1:todel(2)-1, :);
                  N.VERT(todel(2)+1:end, :);
                  NewVert];
        N.TRIV = N.TRIV(N.TRIV(:, 1) ~= N.TRIV(:, 2), :);
        N.TRIV = N.TRIV(N.TRIV(:, 1) ~= N.TRIV(:, 3), :);
        N.TRIV = N.TRIV(N.TRIV(:, 2) ~= N.TRIV(:, 3), :);
        N.n = size(N.VERT, 1);
        N.m = size(N.TRIV, 1);
    end
    
    if isfield(M, 'X')
        N.X = N.VERT(:, 1);
    end
    if isfield(M, 'Y')
        N.Y = N.VERT(:, 2);
    end
    if isfield(M, 'Z')
        N.Z = N.VERT(:, 3);
    end
end

