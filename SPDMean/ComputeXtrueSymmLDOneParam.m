function Xtrue = ComputeXtrueSymmLDOneParam(As, alpha, TOL)


n = size(As{1}, 1);
k = length(As);
Ls = zeros(n, n, k);
for i = 1 : k
    Ls(:, :, i) = As{i};
end


fprintf('I am trying hard to compute the true mean \n');
fprintf('==================================\n');
fprintf('===================\n');
fprintf('=========\n');
%% set general parameters for C++
SolverParams.Tolerance = TOL;
SolverParams.LS_alpha = 1e-4;
SolverParams.LS_beta = 0.999;
SolverParams.LS_ratio1 = 0.5;
SolverParams.LS_ratio2 = 0.5;
SolverParams.Accuracy = 1e-4;
SolverParams.Finalstepsize = -1;
SolverParams.Max_Iteration = 200;
SolverParams.Maxstepsize = 100;
SolverParams.method = 'LRBFGS';
SolverParams.LineSearch_LS = 0; % 'ARMIJO'
SolverParams.Initstepsize = 1;
SolverParams.InitSteptype = 0; % one
SolverParams.HasHHR = 1;
SolverParams.DEBUG = 3;
SolverParams.LengthSY = 2;
SolverParams.isconvex = 0;
SolverParams.IsCheckGradHess = 0;




HasHHR = SolverParams.HasHHR;

Xinitial = computeAH(As, 1);
[Xopt, ~] = TestSPDMeanSymmLDOneParam(Ls, alpha, Xinitial, HasHHR, SolverParams);

Xtrue = reshape(Xopt.main, n, n);

% 
% figure(300)
% D = eig(Xtrue);
% for i = 1 : n
%     refline(0, D(i));
%     hold on;
% end

end