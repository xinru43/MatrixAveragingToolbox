function Xtrue = ComputeXtrueEuc(As, TOL)


n = size(As{1}, 1);
k = length(As);
Ls = zeros(n, n, k);
for i = 1 : k
    Ls(:, :, i) = chol(As{i}, 'lower');
end
    
    
%% set general parameters for C++
SolverParams.Tolerance = TOL;
SolverParams.LS_alpha = 1e-4;
SolverParams.LS_beta = 0.999;
SolverParams.LS_ratio1 = 0.5;
SolverParams.LS_ratio2 = 0.5;
SolverParams.Accuracy = 1e-3;
SolverParams.Finalstepsize = -1;
SolverParams.Max_Iteration = 200;
SolverParams.Maxstepsize = 100;
SolverParams.method = 'LRBFGS';
SolverParams.LineSearch_LS = 0; % 'ARMIJO'
SolverParams.Initstepsize = 1;
SolverParams.InitSteptype = 0; % one
SolverParams.HasHHR = 1;
SolverParams.DEBUG = 1;
SolverParams.LengthSY = 4;
SolverParams.isconvex = 1;
% SolverParams.IsCheckParams = 0;
% SolverParams.IsCheckGradHess = 0;



HasHHR = SolverParams.HasHHR;

Xinitial = computeAH(As, 0);
[Xopt, ~] = TestSPDEucMean(Ls, Xinitial, HasHHR, SolverParams);
Xtrue = reshape(Xopt.main, n, n);


end