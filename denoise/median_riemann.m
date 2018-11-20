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
SolverParams.Accuracy = 1e-4;
SolverParams.Finalstepsize = -1;
SolverParams.Max_Iteration = 50;
SolverParams.Maxstepsize = 100;


SolverParams.method = 'LRBFGS';
SolverParams.LengthSY = 2;
SolverParams.BBratio = 1;
SolverParams.Num_pre_BB = 0;
SolverParams.LineSearch_LS = 4; 
SolverParams.lambdaLower = 1e-2;
SolverParams.lambdaUpper = 1e2;


SolverParams.Initstepsize = 1;
SolverParams.InitSteptype = 0; % one
SolverParams.HasHHR = 1;
SolverParams.isconvex = 1;
SolverParams.DEBUG = 0;




HasHHR = SolverParams.HasHHR;


% compute initial 
mean_type = 'symmldbreg';
[x1, x1time] = mean_closed(As, n, k, mean_type);
x1 = reshape(x1, n, n);


[Xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, comtime] = TestSPDMeanL1(Ls, x1, HasHHR, SolverParams);
comtime = comtime + x1time;

output = reshape(Xopt.main, n, n);


end