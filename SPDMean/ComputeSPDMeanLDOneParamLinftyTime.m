function [times, ComTime] = ComputeSPDMeanLDOneParamLinftyTime(As, alpha, x1, SolverParams, Xtrue)
    
    n = size(As{1}, 1);
    k = length(As);
    Ls = zeros(n, n, k);
    LDAs = zeros(k, 1);
    for i = 1 : k
        Ls(:, :, i) = As{i};
        LDAs(i) = log(det(As{i}));
    end

    HasHHR = SolverParams.HasHHR;
    
    [Xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funs, grads, times] = TestSPDMeanLDOneParamLinfty(Ls, LDAs, alpha, x1, HasHHR, SolverParams, Xtrue);

end
