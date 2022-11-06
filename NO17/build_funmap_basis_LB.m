function [C, FCTSrc, FCTTrg] = build_funmap_basis_LB(Src, Trg, BSrc, BTrg, KSrc_lb, KTrg_lb, SrcLM, TrgLM, para)
    %% default parameters
    % Infer landmarks, if not given
    if nargin < 7
        SrcLM = (500:1000:6000)';   %6 landmarks
        TrgLM = SrcLM;
    end
    LM = [SrcLM, TrgLM];
    
    if nargin >= 9 && isfield(para, 'icp_passes')
        icp_passes = para.icp_passes;
    else
        icp_passes = 5;
    end
    
    if nargin >= 9 && isfield(para, 'maxIter')
        maxIter = para.maxIter;
    else
        maxIter = 500;
    end
    
    if nargin >= 9 && isfield(para, 'a')
        a = para.a;
    else
        a = 1e-1;
    end
    
    if nargin >= 9 && isfield(para, 'b')
        b = para.b;
    else
        b = 1;
    end
    
    if nargin >= 9 && isfield(para, 'c')
        c = para.c;
    else
        c = 1e-3;
    end
        
    
    A_src = Src.A;
    A_trg = Trg.A;
    
    KSrc = size(BSrc, 2);
    KTrg = size(BTrg, 2);
    
    %% Compute descriptors
    FCTSrc = [];
    FCTSrc = [FCTSrc, waveKernelSignature(Src.Phi(:, 1:KSrc_lb), Src.Lambda(1:KSrc_lb), Src.A, 200)];   %200 constraints
    FCTSrc = [FCTSrc, heatKernelMap(Src.Phi(:, 1:KSrc_lb), Src.Lambda(1:KSrc_lb), Src.A, 200, LM(:, 1))];   %1200 constraints
    FCTSrc = FCTSrc(:,1:10:end);

    FCTTrg = [];
    FCTTrg = [FCTTrg, waveKernelSignature(Trg.Phi(:, 1:KTrg_lb), Trg.Lambda(1:KTrg_lb), Trg.A, 200)];
    FCTTrg = [FCTTrg, heatKernelMap(Trg.Phi(:, 1:KTrg_lb), Trg.Lambda(1:KTrg_lb), Trg.A, 200, LM(:, 2))];
    FCTTrg = FCTTrg(:,1:10:end);
    
    %% Normalization
    no = sqrt(diag(FCTSrc' * Src.A * FCTSrc))';
    FCTSrc = FCTSrc ./ repmat(no, [Src.n, 1]);   
    FCTTrg = FCTTrg ./ repmat(no, [Trg.n, 1]);
        
    
    %% Multiplication Operators
    numFct = size(FCTSrc,2);
    OpSrc = cell(numFct, 1);
    OpTrg = cell(numFct, 1);
    for i = 1:numFct
        OpSrc{i} = BSrc' * A_src * (repmat(FCTSrc(:,i), [1, KSrc]) .* BSrc);
        OpTrg{i} = BTrg' * A_trg * (repmat(FCTTrg(:,i), [1, KTrg]) .* BTrg);
    end
    
    Fct_src = BSrc' * A_src * FCTSrc;
    Fct_tar = BTrg' * A_trg * FCTTrg;
    
    %% Laplace operator 
    LbSrc = BSrc' * Src.S * BSrc;   %this comes from the cancellation in BSrc' *A_src * 1./A_src * Src.S * BSrc;
    LbSrc = LbSrc';
    LbTrg = BTrg' * Trg.S * BTrg;   
    LbTrg = LbTrg';
    
    constFct = sign(Src.Phi(1,1) * Trg.Phi(1,1)) * [sqrt(sum(mesh.proc.tri_areas(Trg)) / sum(mesh.proc.tri_areas(Src))); zeros(KTrg - 1, 1)];
    
    %% Fmap Computation
    a = a;      % Descriptors preservation
    b = b;      % Commutativity with descriptors
    c = c;      % Commutativity with Laplacian
    
    funObj = @(F) deal( a * sum(sum((reshape(F, [KTrg, KSrc]) * Fct_src - Fct_tar).^2)) /2 + b * sum(cell2mat(cellfun(@(X,Y) sum(sum((X*reshape(F, [KTrg,KSrc]) - reshape(F, [KTrg,KSrc])*Y).^2)), OpTrg', OpSrc', 'UniformOutput', false)), 2)/2 + c*sum(sum((LbTrg*reshape(F, [KTrg,KSrc]) - reshape(F, [KTrg,KSrc])*LbSrc).^2)),...
                a*vec((reshape(F, [KTrg,KSrc])*Fct_src - Fct_tar)*Fct_src') + b*sum(cell2mat(cellfun(@(X,Y) vec(X'*(X*reshape(F, [KTrg,KSrc]) - reshape(F, [KTrg,KSrc])*Y) - (X*reshape(F, [KTrg,KSrc]) - reshape(F, [KTrg,KSrc])*Y)*Y'), OpTrg', OpSrc', 'UniformOutput', false)), 2) + c*vec(LbTrg'*(LbTrg*reshape(F, [KTrg,KSrc]) - reshape(F, [KTrg,KSrc])*LbSrc) - (LbTrg*reshape(F, [KTrg,KSrc]) - reshape(F, [KTrg,KSrc])*LbSrc)*LbSrc'));
    funProj = @(F) [constFct; F(KTrg+1:end)];
    
    F_lb = zeros(KTrg*KSrc, 1); F_lb(1) = constFct(1);
        
    % Compute the optional functional map using a quasi-Newton method.
    options.verbose = 1;
    options.maxIter = maxIter;
    F_lb = reshape(minConf_PQN(funObj, F_lb, funProj, options), [KTrg,KSrc]);

    %
    if icp_passes
        [C, ~] = icp_refine(BSrc, BTrg, F_lb, icp_passes);
    else
        C = F_lb;
    end

end