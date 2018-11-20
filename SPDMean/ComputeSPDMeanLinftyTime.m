function [times, ComTime] = ComputeSPDMeanLinftyTime(As, x1, SolverParams, Xtrue)
    
    n = size(As{1}, 1);
    k = length(As);
    Ls = zeros(n, n, k);
    for i = 1 : k
        Ls(:, :, i) = chol(As{i}, 'lower');
    end

    HasHHR = SolverParams.HasHHR;
    
    [Xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funs, grads, times] = TestSPDMeanLinfty(Ls, x1, HasHHR, SolverParams, Xtrue);

end
