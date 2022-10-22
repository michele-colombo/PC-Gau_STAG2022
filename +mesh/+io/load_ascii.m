function M = load_ascii(filename)
    % filename end with .ascii
    
    [tempDir, tempFile] = fileparts(filename); 
    tri_file_path = fullfile(tempDir, [tempFile, '.tri']);
    vert_file_path = fullfile(tempDir, [tempFile, '.vert']);
    
    % load tri file
    tri_file = fopen(tri_file_path, 'r');
    A = fscanf(tri_file, '%d %d %d\n', [3 Inf]);
    M.TRIV = A';

    % loade vert file
    vert_file = fopen(vert_file_path, 'r');
    A = fscanf(vert_file, '%f %f %f\n', [3 Inf]);
    M.VERT = A';
    M.X = M.VERT(:,1);
    M.Y = M.VERT(:,2);
    M.Z = M.VERT(:,3);
    
    % number of vertices and faces
    M.n = size(M.VERT, 1);
    M.m = size(M.TRIV, 1);
    
end