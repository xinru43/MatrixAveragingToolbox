function [output, comtime] = median_riemann(As, TOL)


n = size(As(:, :, 1), 1);
k = size(As, 3);
Ls = zeros(n, n, k);
for i = 1 : k
    Ls(:, :, i) = chol(As(:, :, i), 'lower');
end
    
    
%% set general parameters for C++
SolverParams.Tolerance = TOL;
SolverParams.LS_alpha = 1e-4;
SolverParams.LS_beta = 0.999;
SolverParams.LS_ratio1 = 0.5;
SolverParams.LS_ratio2 = 0.5;
SolverParams.Accuracy = 1e-3;
SolverParams.Finalstepsize = -1;
SolverParams.Max_Iteration = 50;
SolverParams.Maxstepsize = 100;
SolverParams.method = 'LRBFGS';
SolverParams.LineSearch_LS = 0; % 'ARMIJO'
SolverParams.Initstepsize = 1;
SolverParams.InitSteptype = 0; % one
SolverParams.HasHHR = 1;
SolverParams.DEBUG = 1;
SolverParams.LengthSY = 2;
SolverParams.isconvex = 0;
HasHHR = SolverParams.HasHHR;


x1 = As(:, :, 1);
for i = 2 : k
    x1 = x1 + As(:, :, i);
end
x1 = x1/k;

[Xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, comtime] = TestSPDMeanL1(Ls, x1, HasHHR, SolverParams);
output = reshape(Xopt.main, n, n);


end