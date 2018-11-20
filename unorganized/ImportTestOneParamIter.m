function ImportTestOneParamIter(As)

% investigate the relationship between the BB stepsize and eigenvalues of
% the Hessian

%% generate data: 100 sets with n = 30, K = 30

n = size(As{1}, 2);
k = length(As);


Ls = zeros(n, n, k);
for i = 1 : k
    Ls(:, :, i) = As{i};
end


%% initial iterate
x1 = computeAH(As, 0); % Type = 0: arithmetic mean; Type = 1: arithmetic-harmonic mean



%% compute Xtrue
alpha = 0.0;

%% set parameters
SolverParams.IsCheckParams = 0;
SolverParams.Stop_Criterion = 2;
SolverParams.Tolerance = 1e-6;
SolverParams.Accuracy = 1e-5;
SolverParams.LS_ratio1 = 0.5;
SolverParams.LS_ratio2 = 0.5;
SolverParams.Minstepsize = eps;
SolverParams.Maxstepsize = 100;
SolverParams.Initstepsize = 1;
SolverParams.Finalstepsize = -1;
SolverParams.LineSearch_LS = 0; % 'ARMIJO'
SolverParams.DEBUG = 0;
SolverParams.HasHHR = 1;
SolverParams.isconvex = 1;
SolverParams.Max_Iteration = 100;

TestMethods = {'BB1', 'BB2', 'ABB', 'FP', 'RNewton'};
BB1 = {'m0WL', 'm0', 'm4', 'm16'};
BB2 = {'m0WL', 'm0', 'm4', 'm16'};
ABB = {'m0WL', 'm0', 'm4', 'm16'};

%% fixed point iteration
if (ismember('FP', TestMethods))
    AsFP = [];
    for i = 1 : k
        AsFP = [AsFP reshape(As{i}, 1, n * n)];
    end
    x1FP = reshape(x1, 1, n * n);
    XtrueFP = reshape(eye(n), 1, n * n);
    Max_Iteration = 200;
    Debug = 2;
    [TFP, ComTimeFP, iterFP, X, disFP, truedisFP, gradsFP, gfgf0FP] = TestSPDDLOneParamFP(AsFP, alpha, x1FP, XtrueFP, n, k, Max_Iteration, 1e-7, Debug);
end
 
%% BB1 family
if (ismember('BB1', TestMethods))
    
    % BB1 parameter
    SolverParams.BBratio = 0;
    SolverParams.Num_pre_BB = 0;
    
    % BB1 without locking condition
    if(ismember('m0WL', BB1))
        SolverParams.method = 'RSD';
        SolverParams.InitSteptype = 1; % BB stepsize
        [gfBB1WL, gfgf0BB1WL, iterBB1WL, nfBB1WL, ngBB1WL, nRBB1WL, nVBB1WL, nVpBB1WL, nHBB1WL, ComTimeBB1WL] = ComputeSPDMeanDLOneParamIter(As, alpha, x1, SolverParams);
    end
    
    
    if(ismember('m0', BB1))
        SolverParams.method = 'LRBFGS';
        SolverParams.LengthSY = 0;
        SolverParams.InitSteptype = 0; % one
        [gfBB1m0, gfgf0BB1m0, iterBB1m0, nfBB1m0, ngBB1m0, nRBB1m0, nVBB1m0, nVpBB1m0, nHBB1m0, ComTimeBB1m0] = ComputeSPDMeanDLOneParamIter(As, alpha, x1, SolverParams);
    end
    
    
    if(ismember('m2', BB1))
        SolverParams.method = 'LRBFGS';
        SolverParams.LengthSY = 2;
        SolverParams.InitSteptype = 0; % one
        [gfBB1m4, gfgf0BB1m4, iterBB1m4, nfBB1m4, ngBB1m4, nRBB1m4, nVBB1m4, nVpBB1m4, nHBB1m4, ComTimeBB1m4] = ComputeSPDMeanDLOneParamIter(As, alpha, x1, SolverParams);
    end
    
    
    if(ismember('m4', BB1))
        SolverParams.method = 'LRBFGS';
        SolverParams.LengthSY = 4;
        SolverParams.InitSteptype = 0; % one
        [gfBB1m4, gfgf0BB1m4, iterBB1m4, nfBB1m4, ngBB1m4, nRBB1m4, nVBB1m4, nVpBB1m4, nHBB1m4, ComTimeBB1m4] = ComputeSPDMeanDLOneParamIter(As, alpha, x1, SolverParams);
    end
    
    
    if(ismember('m16', BB1))
        SolverParams.method = 'LRBFGS';
        SolverParams.LengthSY = 16;
        SolverParams.InitSteptype = 0; % one
        [gfBB1m16, gfgf0BB1m16, iterBB1m16, nfBB1m16, ngBB1m16, nRBB1m16, nVBB1m16, nVpBB1m16, nHBB1m16, ComTimeBB1m16] = ComputeSPDMeanDLOneParamIter(As, alpha, x1, SolverParams);
    end
    
end



%% BB2 family
if (ismember('BB2', TestMethods))
    
    % BB2 parameter
    SolverParams.BBratio = 1;
    SolverParams.Num_pre_BB = 0;
    
    % BB2 without locking condition
    if(ismember('m0WL', BB2))
        SolverParams.method = 'RSD';
        SolverParams.InitSteptype = 1; % BB stepsize
        [gfBB2WL, gfgf0BB2WL, iterBB2WL, nfBB2WL, ngBB2WL, nRBB2WL, nVBB2WL, nVpBB2WL, nHBB2WL, ComTimeBB2WL] = ComputeSPDMeanDLOneParamIter(As, alpha, x1, SolverParams);
    end
    
    
    if(ismember('m0', BB2))
        SolverParams.method = 'LRBFGS';
        SolverParams.LengthSY = 0;
        SolverParams.InitSteptype = 0; % one
        [gfBB2m0, gfgf0BB2m0, iterBB2m0, nfBB2m0, ngBB2m0, nRBB2m0, nVBB2m0, nVpBB2m0, nHBB2m0, ComTimeBB2m0] = ComputeSPDMeanDLOneParamIter(As, alpha, x1, SolverParams);
    end
    
    
    if(ismember('m2', BB2))
        SolverParams.method = 'LRBFGS';
        SolverParams.LengthSY = 2;
        SolverParams.InitSteptype = 0; % one
        [gfBB1m4, gfgf0BB1m4, iterBB1m4, nfBB1m4, ngBB1m4, nRBB1m4, nVBB1m4, nVpBB1m4, nHBB1m4, ComTimeBB1m4] = ComputeSPDMeanDLOneParamIter(As, alpha, x1, SolverParams);
    end
    
    
    if(ismember('m4', BB2))
        SolverParams.method = 'LRBFGS';
        SolverParams.LengthSY = 4;
        SolverParams.InitSteptype = 0; % one
        [gfBB2m4, gfgf0BB2m4, iterBB2m4, nfBB2m4, ngBB2m4, nRBB2m4, nVBB2m4, nVpBB2m4, nHBB2m4, ComTimeBB2m4] = ComputeSPDMeanDLOneParamIter(As, alpha, x1, SolverParams);
    end
    
    
    if(ismember('m16', BB2))
        SolverParams.method = 'LRBFGS';
        SolverParams.LengthSY = 16;
        SolverParams.InitSteptype = 0; % one
        [gfBB2m16, gfgf0BB2m16, iterBB2m16, nfBB2m16, ngBB2m16, nRBB2m16, nVBB2m16, nVpBB2m16, nHBB2m16, ComTimeBB2m16] = ComputeSPDMeanDLOneParamIter(As, alpha, x1, SolverParams);
    end
    
end






%% ABB family
if (ismember('ABB', TestMethods))
    
    % ABB parameter
    SolverParams.BBratio = 0.5;
    SolverParams.Num_pre_BB = 10;
    
    % ABB without locking condition
    if(ismember('m0WL', ABB))
        SolverParams.method = 'RSD';
        SolverParams.InitSteptype = 1; % BB stepsize
        [gfABBWL, gfgf0ABBWL, iterABBWL, nfABBWL, ngABBWL, nRABBWL, nVABBWL, nVpABBWL, nHABBWL, ComTimeABBWL] = ComputeSPDMeanDLOneParamIter(As, alpha, x1, SolverParams);
    end
    
    
    if(ismember('m0', ABB))
        SolverParams.method = 'LRBFGS';
        SolverParams.LengthSY = 0;
        SolverParams.InitSteptype = 0; % one
        [gfABBm0, gfgf0ABBm0, iterABBm0, nfABBm0, ngABBm0, nRABBm0, nVABBm0, nVpABBm0, nHABBm0, ComTimeABBm0] = ComputeSPDMeanDLOneParamIter(As, alpha, x1, SolverParams);
    end
    
    
   if(ismember('m2', ABB))
        SolverParams.method = 'LRBFGS';
        SolverParams.LengthSY = 2;
        SolverParams.InitSteptype = 0; % one
        [gfBB1m4, gfgf0BB1m4, iterBB1m4, nfBB1m4, ngBB1m4, nRBB1m4, nVBB1m4, nVpBB1m4, nHBB1m4, ComTimeBB1m4] = ComputeSPDMeanDLOneParamIter(As, alpha, x1, SolverParams);
    end
    
    
    if(ismember('m4', ABB))
        SolverParams.method = 'LRBFGS';
        SolverParams.LengthSY = 4;
        SolverParams.InitSteptype = 0; % one
        [gfABBm4, gfgf0ABBm4, iterABBm4, nfABBm4, ngABBm4, nRABBm4, nVABBm4, nVpABBm4, nHABBm4, ComTimeABBm4] = ComputeSPDMeanDLOneParamIter(As, alpha, x1, SolverParams);
    end
    
    
    if(ismember('m16', ABB))
        SolverParams.method = 'LRBFGS';
        SolverParams.LengthSY = 16;
        SolverParams.InitSteptype = 0; % one
        [gfABBm16, gfgf0ABBm16, iterABBm16, nfABBm16, ngABBm16, nRABBm16, nVABBm16, nVpABBm16, nHABBm16, ComTimeABBm16] = ComputeSPDMeanDLOneParamIter(As, alpha, x1, SolverParams);
    end
    
end



    if(ismember('RNewton', TestMethods))
        SolverParams.method = 'RNewton';
        SolverParams.InitSteptype = 0; % one
        [gfRNewton, gfgf0RNewton, iterRNewton, nfRNewton, ngRNewton, nRRNewton, nVRNewton, nVpRNewton, nHRNewton, ComTimeRNewton] = ComputeSPDMeanDLOneParamIter(As, alpha, x1, SolverParams);
    end

end




