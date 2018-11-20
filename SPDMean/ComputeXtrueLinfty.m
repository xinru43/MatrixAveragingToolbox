function Xtrue = ComputeXtrueLinfty(As, Xinitial, TOL)


n = size(As{1}, 1);
k = length(As);
Ls = zeros(n, n, k);
for i = 1 : k
    Ls(:, :, i) = chol(As{i}, 'lower');
end


%% set general parameters for C++
SolverParams.Tolerance = TOL;
SolverParams.Accuracy = eps;
SolverParams.Finalstepsize = 1;
SolverParams.Stop_Criterion = 3;
SolverParams.LS_alpha = 1e-4;
SolverParams.LS_beta = 0.999;
SolverParams.Minstepsize = eps;
SolverParams.Maxstepsize = 200;
SolverParams.Initstepsize = 1;


SolverParams.lambdaLower = 1e-2;
SolverParams.lambdaUpper = 50;
SolverParams.Eps = 1e-4;
SolverParams.Del = 1e-6;
SolverParams.Theta_del = 1e-2;
SolverParams.Theta_eps = 1e-2;


SolverParams.Max_Iteration = 200;
SolverParams.LineSearch_LS = 4; 
SolverParams.InitSteptype = 0; 
SolverParams.isconvex = 1;
SolverParams.HasHHR = 1;
SolverParams.DEBUG = 1;


SolverParams.method = 'LRBFGSLPSub';
SolverParams.LengthSY = 4;



HasHHR = 1;
[Xopt, ~] = TestSPDMeanLinfty(Ls, Xinitial, HasHHR, SolverParams);
Xtrue = reshape(Xopt.main, n, n);


end