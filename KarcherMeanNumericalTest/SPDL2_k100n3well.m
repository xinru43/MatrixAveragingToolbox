function SPDL2_k100n3well

% number of data tested
Ndata = 50;

% specify the path of data and results
DataPath = sprintf('L2Results/k100n3well/');
SavePath = sprintf('L2Results/k100n3well/');

% methods tested
TestMethods = {'RSDQR', 'LRBFGS', 'BINI', 'RBFGS'};
metric = 'Rie';

% memory size for LRBFGS
MemorySize = {0, 2, 4};

% stabletime * Ndata = 1 minute
stabletime = 1.2;
stabletime = 1e-4;

% used for average distance
dt = 1e-5;


% =============================================================================
% Run Experiments
% =============================================================================
rng(1);
TestSeed = 10000 * rand(Ndata, 1);
save([SavePath 'TestSeed.mat'], 'TestSeed');

for times = 1 : Ndata
    
    fprintf('======================== data: %d\n', times);
    
    %% start the timer to remove the overhead
    TestTimerStart;
    
    
    %% import data
    generate_new_data = 0;
    
    if generate_new_data
        n = 3;
        k = 100;
        f = 1; % controls the condition number of data points
        [Kmean, As, cond_Kmean, cond_As] = generate_mean_data(k, n, f, TestSeed(times));
        save([SavePath 'Kmean' num2str(times) '.mat'], 'As', 'Kmean', 'cond_As', 'cond_Kmean');
        Xtrue = Kmean;
    else
        Data = importdata([DataPath 'Kmean' num2str(times) '.mat']);
        As = Data.As;
        n = size(As{1}, 2);
        k = length(As);
        Xtrue = Data.Kmean;
    end
    
    
    % Compute the true solution when necessary
    if ~exist('Xtrue')
        TOL = 1e-16;
        Xtrue = ComputeXtrue(As, x1, TOL);
    end
    

    % compute the initial point
    x1 = computeAH(As, 1);
    
    
    %% never change parameters
    SolverParams.IsCheckParams = 0;
    SolverParams.IsCheckGradHess = 0;
    SolverParams.LS_alpha = 1e-4;
    SolverParams.LS_beta = 0.999;
    SolverParams.Stop_Criterion = 3;
    SolverParams.Tolerance = 1e-17;
    SolverParams.Minstepsize = eps;
    SolverParams.Maxstepsize = 10;
    SolverParams.Max_Iteration = 50;
    SolverParams.isconvex = 1;
    SolverParams.LineSearch_LS = 0; % Armijo 
    
    
    SolverParams.Accuracy = 1e-5;
    SolverParams.HasHHR = 1;
    SolverParams.LS_ratio1 = 0.5;
    SolverParams.LS_ratio2 = 0.5;    
    SolverParams.Finalstepsize = -1;
    
    
    % compute initial step size in the first iteration: h0 = 2/(1 + B)
    c = zeros(1, k);
    L = chol(x1, 'lower');
    for i = 1 : k
        c(i) = cond((L\As{i})/L');
    end
    QRstepsize = sum(((c + 1).*log(c))./(c - 1));
    QRstepsize = QRstepsize/(2.0 * k);
    QRstepsize = 2.0/(1.0 + QRstepsize);
    SolverParams.Initstepsize = QRstepsize;
    
    
     %% =========================== RBB ===========================
    if (ismember('RSDBB', TestMethods))
        SolverParams.method = 'RSD';
        SolverParams.InitSteptype = 1; % use BB s's/s'y version
        SolverParams.BBratio = 1;
        SolverParams.Num_pre_BB = 0;
        inFileName = 'RSDBB';
        TestPipeL2(As, x1, Xtrue, SolverParams, stabletime, SavePath, inFileName, times, metric)
    end
    
    
    %% =========================== LRBFGS-BB2 ===========================
    if (ismember('LRBFGS', TestMethods))
        for m = 1 : length(MemorySize)
            SolverParams.method = 'LRBFGS';
            SolverParams.InitSteptype = 0; % one
            SolverParams.BBratio = 1;
            SolverParams.Num_pre_BB = 0;
            SolverParams.LengthSY = MemorySize{m};            
            inFileName = ['LRBFGSm' int2str(MemorySize{m})];
            TestPipeL2(As, x1, Xtrue, SolverParams, stabletime, SavePath, inFileName, times, metric)
        end
    end
    
        
    %% =========================== RBFGS ===========================
    if (ismember('RBFGS', TestMethods))
        SolverParams.method = 'RBFGS';
        SolverParams.InitSteptype = 0; % one
        inFileName = 'RBFGS';
        TestPipeL2(As, x1, Xtrue, SolverParams, stabletime, SavePath, inFileName, times, metric)
    end 
       
    
    
    %% =========================== RSD + QR stepsize + without linesearch ===========================
    if (ismember('RSDQR', TestMethods))
        inFileName = 'RSDQR';
        Max_Iteration = 50;
        TOL = SolverParams.Tolerance;
        Retraction = 2;
        StepType = 0;
        TestPipeL2RSD(As, x1, Xtrue, Max_Iteration, TOL, Retraction, StepType, stabletime, SavePath, inFileName, times)
    end
    
   
     %% =========================== Bini ===========================
    if (ismember('BINI', TestMethods))
        inFileName = 'BINI';
        Max_Iteration = 50;
        TOL = SolverParams.Tolerance;
        
        Retraction = 1;
        StepType = 1;
        TestPipeL2RSD(As, x1, Xtrue, Max_Iteration, TOL, Retraction, StepType, stabletime, SavePath, inFileName, times)
    end
    
end




%% =============================================================================
% Average Data
% =============================================================================

if (ismember('RBFGS', TestMethods))
    inFileName = 'RBFGS';
    FileName = ['Ave' inFileName];
    for times = 1 : Ndata
        inData{times} = importdata([SavePath inFileName int2str(times) '.mat']);
        delete([SavePath inFileName num2str(times) '.mat']);
    end
     AveDataL2(inData, dt, SavePath, FileName);
end


if (ismember('RSDBB', TestMethods))
    inFileName = 'RSDBB';
    FileName = ['Ave' inFileName];
    for times = 1 : Ndata
        inData{times} = importdata([SavePath inFileName int2str(times) '.mat']);
        delete([SavePath inFileName num2str(times) '.mat']);
    end
     AveDataL2(inData, dt, SavePath, FileName);
end


if (ismember('LRBFGS', TestMethods))
    for m = 1 : length(MemorySize)
        inFileName = ['LRBFGSm' int2str(MemorySize{m})];
        FileName = ['Ave' inFileName];
        
        for times = 1 : Ndata
            inData{times} = importdata([SavePath inFileName int2str(times) '.mat']);
            delete([SavePath inFileName num2str(times) '.mat']);
        end
        AveDataL2(inData, dt, SavePath, FileName);
     end
end


if (ismember('RSDQR', TestMethods))
    inFileName = 'RSDQR';
    FileName = ['Ave' inFileName];
    for times = 1 : Ndata
        inData{times} = importdata([SavePath inFileName int2str(times) '.mat']);
        delete([SavePath inFileName num2str(times) '.mat']);
    end
    AveDataRSD(inData, dt, SavePath, FileName);
end



if (ismember('BINI', TestMethods))
    inFileName = 'BINI';
    FileName = ['Ave' inFileName];
    for times = 1 : Ndata
        inData{times} = importdata([SavePath inFileName int2str(times) '.mat']);
        delete([SavePath inFileName num2str(times) '.mat']);
    end
    AveDataRSD(inData, dt, SavePath, FileName);
end

end







