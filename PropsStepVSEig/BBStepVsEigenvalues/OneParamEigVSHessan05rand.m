function OneParamEigVSHessan05rand


Ndata = 1;


TestMethods = {'BB1', 'BB2', 'ABB', 'LRBFGS', 'RSD'};

stabletime = 1000;
alpha = -0.5;
dt = 1e-5;

AsPath = sprintf('PropsStepVSEig/');
DataPath = sprintf('PropsStepVSEig/an05rand/');
SavePath = sprintf('/panfs/storage.local/home-4/xy12d/ToXinru/PropsStepVSEig/an05rand/');


RSDBBiter = 60;
RSDABBiter = 60;
LRBFGSiter = 60;
RBFGSiter = 60;
RNewtoniter = 60;
RSDiter = 60;


% =============================================================================
% Run Experiments
% =============================================================================

for times = 1 : Ndata
    
    %% start the timer to remove the overhead
    TestTimerStart;
    
    %% import data including initial point
    As = importdata([AsPath 'As.mat']);
    n = size(As{1}, 2);
    k = length(As);
    
    
    x1 = importdata([AsPath 'x1an05.mat']);
    
    
    
    %% never change parameters
SolverParams.IsCheckParams = 0;
SolverParams.Tolerance = 1e-17;
SolverParams.Accuracy = 1e-4;
SolverParams.LS_ratio1 = 0.5;
SolverParams.LS_ratio2 = 0.5;
SolverParams.Minstepsize = eps;
SolverParams.Maxstepsize = 100;
SolverParams.Initstepsize = 1;
SolverParams.Finalstepsize = -1;
SolverParams.LineSearch_LS = 0; % 'ARMIJO'
SolverParams.isconvex = 0;
SolverParams.Max_Iteration = 60;
SolverParams.HasHHR = 1;
    
    
    
    
    %% =========================== Compute the true solution ===========================
    Xtrue = ComputeXtrueDLOneParam(As, alpha, 1e-12);
    
    
    
    %% =========================== BB1 ===========================
    if (ismember('BB1', TestMethods))
                SolverParams.BBratio = 0;
                SolverParams.Num_pre_BB = 0;
                SolverParams.method = 'RSD';
                SolverParams.InitSteptype = 1;
                inFileName = ['BB1'];
                TestPipe(As, alpha, x1, Xtrue, SolverParams, stabletime, SavePath, inFileName, times)
    end
    
    
   %% =========================== BB2 ===========================
    if (ismember('BB2', TestMethods))
                SolverParams.BBratio = 1;
                SolverParams.Num_pre_BB = 0;
                SolverParams.method = 'RSD';
                SolverParams.InitSteptype = 1;
                inFileName = ['BB2'];
                TestPipe(As, alpha, x1, Xtrue, SolverParams, stabletime, SavePath, inFileName, times)
    end
    
    
    %% =========================== ABB ===========================
    if (ismember('ABB', TestMethods))
                SolverParams.BBratio = 0.8;
                SolverParams.Num_pre_BB = 5;
                SolverParams.method = 'RSD';
                SolverParams.InitSteptype = 1;
                inFileName = ['ABB'];
                TestPipe(As, alpha, x1, Xtrue, SolverParams, stabletime, SavePath, inFileName, times)
    end

    
    %% =========================== LRBFGS-BB2 ===========================
    
    if (ismember('LRBFGS', TestMethods))
            SolverParams.BBratio = 1;
            SolverParams.Num_pre_BB = 0;
            SolverParams.LengthSY = 2;
            SolverParams.method = 'LRBFGS';   
            SolverParams.InitSteptype = 0;
            inFileName = 'LRBFGS';
            TestPipe(As, alpha, x1, Xtrue, SolverParams, stabletime, SavePath, inFileName, times)
    end
    

    
    %% =========================== RSD ===========================
    if (ismember('RSD', TestMethods))
        SolverParams.method = 'RSD';
        SolverParams.InitSteptype = 3; 
        inFileName = 'RSD';
        TestPipe(As, alpha, x1, Xtrue, SolverParams, stabletime, SavePath, inFileName, times)
    end
    

end




%% =============================================================================
% Average Data
% =============================================================================

if (ismember('BB1', TestMethods))
    inFileName = 'BB1';
    FileName = ['Ave' inFileName];
    for times = 1 : Ndata
        inData{times} = importdata([SavePath inFileName int2str(times) '.mat']);
        delete([SavePath inFileName int2str(times) '.mat']);
    end
    AveData(inData, dt, stabletime, SavePath, FileName);
end

if (ismember('BB2', TestMethods))
    inFileName = 'BB2';
    FileName = ['Ave' inFileName];
    for times = 1 : Ndata
        inData{times} = importdata([SavePath inFileName int2str(times) '.mat']);
        delete([SavePath inFileName int2str(times) '.mat']);
    end
    AveData(inData, dt, stabletime, SavePath, FileName);
end

if (ismember('ABB', TestMethods))
    inFileName = 'ABB';
    FileName = ['Ave' inFileName];
    for times = 1 : Ndata
        inData{times} = importdata([SavePath inFileName int2str(times) '.mat']);
        delete([SavePath inFileName int2str(times) '.mat']);
    end
    AveData(inData, dt, stabletime, SavePath, FileName);
end


if (ismember('LRBFGS', TestMethods))
    inFileName = 'LRBFGS';
    FileName = ['Ave' inFileName];
    for times = 1 : Ndata
        inData{times} = importdata([SavePath inFileName int2str(times) '.mat']);
        delete([SavePath inFileName int2str(times) '.mat']);
    end
    AveData(inData, dt, stabletime, SavePath, FileName);
end


if (ismember('RSD', TestMethods))
    inFileName = 'RSD';
    FileName = ['Ave' inFileName];
    for times = 1 : Ndata
        inData{times} = importdata([SavePath inFileName int2str(times) '.mat']);
        delete([SavePath inFileName int2str(times) '.mat']);
    end
    AveData(inData, dt, stabletime, SavePath, FileName);
end


end







