function [output, comtime] = minimaxcenter_alpha(As, alpha, TOL)
%     [diss, eucdiss, initstepsize, stepsize, grads, hesseig, Eig, times, X, funs, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime]


n = size(As(:, :, 1), 1);
k = size(As, 3);
   

SolverParams.Tolerance = TOL;
SolverParams.LS_alpha = 1e-4;
SolverParams.LS_beta = 0.999;
SolverParams.LS_ratio1 = 0.5;
SolverParams.LS_ratio2 = 0.5;
SolverParams.Accuracy = 1e-4;
SolverParams.Finalstepsize = -1;
SolverParams.Max_Iteration = 50;
SolverParams.Maxstepsize = 100;


SolverParams.lambdaLower = 1e-2;
SolverParams.lambdaUpper = 1e2;
SolverParams.Eps = 1e-4;
SolverParams.Del = 1e-4;
SolverParams.method = 'RBFGS';
SolverParams.LineSearch_LS = 4; 




SolverParams.Initstepsize = 1;
SolverParams.InitSteptype = 0; % one
SolverParams.HasHHR = 1;
SolverParams.DEBUG = 0;
SolverParams.isconvex = 1;

HasHHR = SolverParams.HasHHR;

% compute initial 
mean_type = 'symmldbreg';
[x1, x1time] = mean_closed(As, n, k, mean_type);
x1 = reshape(x1, n, n);


LDAs = zeros(k, 1);
for i = 1 : k
    LDAs(i) = log(det(As(:, :, i)));
end


% compute the mean
[Xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, comtime] = TestSPDMeanLDOneParamLinfty(As, LDAs, alpha, x1, HasHHR, SolverParams);
comtime = comtime + x1time;



output = reshape(Xopt.main, n, n);
end


