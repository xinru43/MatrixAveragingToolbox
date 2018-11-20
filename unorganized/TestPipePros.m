function TestPipePros(As, alpha, x1, Xtrue, SolverParams, stabletime, SavePath, inFileName, times)

% compute averaged computation time for each data set
SolverParams.DEBUG = 2;
[T0, ComTime0] = ComputeOneParamProsTime(As, alpha, x1, SolverParams, Xtrue);

T = zeros(size(T0));
for st = 1 : stabletime
    [TT, ComTime(st)] = ComputeOneParamProsTime(As, alpha, x1, SolverParams, Xtrue);
    T = T + TT;
end    
T = T/stabletime;
AveComTime = mean(ComTime);

% compute distance
SolverParams.DEBUG = 3;
[~, ~, iter, F, G, dis, initstepsize, stepsize, Eig, f, gf, gfgf0, nf, ng, nR, nV, nVp, nH, X] = ComputeOneParamPros(As, alpha, x1, SolverParams, Xtrue);
  

T = T';
G = G';
dis = dis';
    
save([SavePath inFileName int2str(times) '.mat'],...
    'f', 'gf', 'gfgf0', 'nf', 'ng', 'nR', 'nV', 'nVp', 'nH', 'ComTime', 'AveComTime', ...
    'F','G','T','dis', 'iter', 'initstepsize', 'stepsize', 'T0', 'ComTime0', 'Eig');

end




