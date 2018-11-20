function [times, ComTime, iter, funs, grads, diss, initstepsize, stepsize, Eig, f, gf, gfgf0, nf, ng, nR, nV, nVp, nH, X] = ComputeSPDMeanLDOneParamLinfty(As, alpha, x1, SolverParams, Xtrue)
    
    n = size(As{1}, 1);
    k = length(As);
    Ls = zeros(n, n, k);
    LDAs = zeros(k, 1);
    for i = 1 : k
        Ls(:, :, i) = As{i};
        LDAs(i) = log(det(As{i}));
    end

    HasHHR = SolverParams.HasHHR;
    
    [Xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funs, grads, times, diss, Eig, initstepsize, stepsize] = TestSPDMeanLDOneParamLinfty(Ls, LDAs, alpha, x1, HasHHR, SolverParams, Xtrue);
    X = reshape(Xopt.main, n, n);
end
