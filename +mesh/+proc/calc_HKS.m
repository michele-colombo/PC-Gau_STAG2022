function [HKS,t] = calc_HKS(M, nstep, evecs_size)
    if nargin<3
        evecs_size=size(M.evecs, 2);
    end
    if nargin<2
        nstep = 100;
    end
    nstep = nstep-1;
    evals = M.evals(1:evecs_size);
    evecs = M.evecs(:,1:evecs_size);

    % based on J.Sun's code
    tmin = abs(4*log(10) / evals(end));
    tmax = abs(4*log(10) / evals(2));
    stepsize = (log(tmax) - log(tmin)) / nstep;
    logts = log(tmin):stepsize:log(tmax);
    t = exp(logts);
    %t = 2.^(1: 1/16 : 18-1/16);

    HKS = zeros(M.n, length(t));

    % skipping the first freq. as the authors do in their code
    for i=1:length(t)
        HKS(:,i) = sum(...
            (evecs(:,2:end)).^2 .* repmat(exp(-t(i)*evals(2:end))',M.n,1), ...
            2);
    end

end
