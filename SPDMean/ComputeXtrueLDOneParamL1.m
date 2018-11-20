function Xtrue = ComputeXtrueLDOneParamL1(As, alpha, x1, TOL)


n = size(As{1}, 1);
k = length(As);
Ls = zeros(n, n, k);
LDAs = zeros(k, 1);
for i = 1 : k
    Ls(:, :, i) = As{i};
    LDAs(i) = log(det(As{i}));
end



%% set general parameters for C++
SolverParams.Tolerance = TOL;
SolverParams.Accuracy = 1e-5;
SolverParams.Finalstepsize = 1;
SolverParams.Stop_Criterion = 3;
SolverParams.LS_alpha = 1e-4;
SolverParams.LS_beta = 0.999;
SolverParams.Minstepsize = eps;
SolverParams.Maxstepsize = 200;
SolverParams.Initstepsize = 1;


SolverParams.lambdaLower = 1e-2;
SolverParams.lambdaUpper = 50;
SolverParams.Eps = 1e-2;
SolverParams.Del = 1e-2;
SolverParams.Theta_del = 1e-4;
SolverParams.Theta_eps = 1e-2;


SolverParams.Max_Iteration = 200;
SolverParams.LineSearch_LS = 4; 
SolverParams.InitSteptype = 0; 
SolverParams.isconvex = 1;
SolverParams.HasHHR = 1;
SolverParams.DEBUG = 1;
    

SolverParams.method = 'LRBFGS';
SolverParams.LengthSY = 0;

HasHHR = 1;

[Xopt, ~] = TestSPDMeanLDOneParamL1(Ls, LDAs, alpha, x1, HasHHR, SolverParams);
Xtrue = reshape(Xopt.main, n, n);


end