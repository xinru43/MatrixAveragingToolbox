function output = ComputeOneParamIter(As, alpha, x1, SolverParams, stabletime)
    
    n = size(As{1}, 1);
    k = length(As);
    Ls = zeros(n, n, k);
    for i = 1 : k
        Ls(:, :, i) = As{i};
    end

     
    HasHHR = SolverParams.HasHHR;
    
    % initialize the library
    TestSPDMeanLDOneParam(Ls, alpha, x1, HasHHR, SolverParams);
    
    % take average time
    for st = 1 : stabletime
    [Xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, comtime(st)] = TestSPDMeanLDOneParam(Ls, alpha, x1, HasHHR, SolverParams);
    end
    ComTime = mean(comtime);
    
    output.gf = gf;
    output.gfgf0 = gfgf0;
    output.iter = iter;
    output.nf = nf;
    output.ng = ng;
    output.nR = nR;
    output.nV = nV;
    output.nVp = nVp;
    output.nH = nH;
    output.ComTime = ComTime;
    output.TimePerIter = ComTime/iter;
end