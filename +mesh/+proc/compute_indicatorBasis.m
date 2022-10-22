function Sbasis = compute_indicatorBasis(S,taus,freq,points,thrs, binary, max_eigs)
    if nargin < 7
        max_eigs = length(S.Lambda);
    end
    if nargin < 6
        binary = false;
    end
    if nargin < 5
        thrs = ones(size(taus));
    end
    
    Lambda = S.Lambda(1:max_eigs);
    Phi = S.Phi(:,1:max_eigs);
    
    S_sqrt_area = sqrt(sum(mesh.proc.tri_areas(S)));
    Sbasis = [];

    for t = 1:length(taus)

        % prepare spectral Gaussians
        tau = taus(t)*((S_sqrt_area)^2);
        ghat = exp(-(tau*abs(-Lambda(freq)+Lambda)));
        ghat = ghat./norm(ghat,2);
        coeff = bsxfun(@times,ghat,Phi(points,:)');    

        % compute and normalize spectral Gaussians
        S_basis = Phi * coeff;
        S_basis = bsxfun(@times,abs(S_basis), 1./max(abs(S_basis)));

        % binarization
        if binary 
            thr = thrs(t);
            S_basis(S_basis>thr) = 1;
            S_basis(S_basis<1) = 0;
        end

        Sbasis = [Sbasis, S_basis];

    end

end