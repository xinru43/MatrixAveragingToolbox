function output = ComputeSPDMeanSymmLDOneParam(As, alpha, x1, SolverParams, Xtrue)
%     [diss, eucdiss, initstepsize, stepsize, grads, hesseig, Eig, times, X, funs, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime]
    n = size(As{1}, 1);
    k = length(As);
    Ls = zeros(n, n, k);
    for i = 1 : k
        Ls(:, :, i) = As{i};
    end

     
    HasHHR = SolverParams.HasHHR;
    
    [Xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funs, grads, times, diss, Eig, initstepsize, stepsize, eucdiss, hesseig] = TestSPDMeanSymmLDOneParam(Ls, alpha, x1, HasHHR, SolverParams, Xtrue);
    X = reshape(Xopt.main, n, n);
    output.funs = funs;
    output.grads = grads;
    output.times = times;
    output.diss = diss;
    output.Eig = Eig;
    output.initstepsize = initstepsize;
    output.stepsize = stepsize;
    output.eucdiss = eucdiss;
    output.hesseig = hesseig;
    output.iter = iter;
    output.X = X;
%     output.f = f;
%     output.gf = gf;
%     output.gfgf0 = gfgf0;
%     output.iter = iter;
%     output.nf = nf;
%     output.ng = ng;
%     output.nR = nR;
    
    
end


