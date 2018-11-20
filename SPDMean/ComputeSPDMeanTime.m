function [times, ComTime] = ComputeSPDMeanTime(As, x1, SolverParams, Xtrue, metric)
    
    n = size(As{1}, 1);
    k = length(As);
    Ls = zeros(n, n, k);
    for i = 1 : k
         Ls(:, :, i) = chol(As{i}, 'lower');
    end

    HasHHR = SolverParams.HasHHR;
    
    if (metric == 'Euc')
       [Xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funs, grads, times] = TestSPDEucMean(Ls, x1, HasHHR, SolverParams, Xtrue);
    else
       [Xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funs, grads, times] = TestSPDMean(Ls, x1, HasHHR, SolverParams, Xtrue);
    end



end
