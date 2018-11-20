function [times, ComTime, iter, funs, grads, diss, initstepsize, stepsize, Eig, f, gf, gfgf0, nf, ng, nR, nV, nVp, nH, X] = ComputeSPDMeanL1(As, x1, SolverParams, Xtrue)
    
    n = size(As{1}, 1);
    k = length(As);
    Ls = zeros(n, n, k);
    for i = 1 : k
        Ls(:, :, i) = chol(As{i}, 'lower');
    end

    HasHHR = SolverParams.HasHHR;
    
    [Xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funs, grads, times, diss, Eig, initstepsize, stepsize] = TestSPDMeanL1(Ls, x1, HasHHR, SolverParams, Xtrue);
    X = reshape(Xopt.main, n, n);
end
