function TestPipeL2(As, x1, Xtrue, SolverParams, runingtime, SavePath, inFileName, times, metric)

%% compute averaged computation time for each data set
SolverParams.DEBUG = 2;
[T0, ComTime0] = ComputeSPDMeanTime(As, x1, SolverParams, Xtrue, metric);

T = zeros(size(T0));
ComTime(1) = 0;
st = 1;

while sum(ComTime) < runingtime
    [TT, ComTime(st)] = ComputeSPDMeanTime(As, x1, SolverParams, Xtrue, metric);
    T = T + TT;
    st = st + 1;
end    
T = T/(st - 1);
AveComTime = mean(ComTime);


%% compute distance
SolverParams.DEBUG = 3;
[~, ~, iter, F, G, dis, initstepsize, stepsize, Eig, f, gf, gfgf0, nf, ng, nR, nV, nVp, nH, X] = ComputeSPDMean(As, x1, SolverParams, Xtrue, metric);
  

T = T';
G = G';
dis = dis';
F = F';
    
save([SavePath inFileName int2str(times) '.mat'],...
    'f', 'gf', 'gfgf0', 'nf', 'ng', 'nR', 'nV', 'nVp', 'nH', 'ComTime', 'AveComTime', ...
    'F','G','T','dis', 'iter', 'initstepsize', 'stepsize', 'T0', 'ComTime0', 'Eig');

end




