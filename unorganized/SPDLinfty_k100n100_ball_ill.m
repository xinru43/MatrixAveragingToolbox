function SPDLinfty_k100n3_cluster

cd ToXinru;
ImportDIR;

Ndata = 50;
rng(1);
TestSeed = 10000 * rand(Ndata, 1);


TestMethods = {'AN', 'LRBFGS-WolfeLP', 'RBFGS-WolfeLP'};



stabletime = 1.2;

MemorySize = {0, 2};
dt = 0.5e-4;



DataPath = sprintf('L1Results/k100n100_ball_ill/');
SavePath = sprintf('/gpfs/home/xy12d/ToXinru/LinftyResults/k100n100_ball_ill/');


LRBFGSiter = 60;
RBFGSiter = 60;
RGSiter = 60;



% =============================================================================
% Run Experiments
% =============================================================================


for times = 1 : Ndata
    
    %% start the timer to remove the overhead
    TestTimerStart;
    
    %% import data including initial point
    As = importdata([DataPath 'As' num2str(times) '.mat']);
    n = size(As{1}, 2);
    k = length(As);
    d = (n * (n + 1))/2;
    
    
    %type = 2;
    %As = LinftyData(k, n, type, TestSeed(times));
    %save([SavePath 'As' num2str(times) '.mat'], 'As');

    
    x1 = As{1};
    
    
    %% never change parameters
    SolverParams.Tolerance = 1e-17;
    SolverParams.Accuracy = eps;
    SolverParams.Finalstepsize = 1;
    SolverParams.Stop_Criterion = 3;
    SolverParams.LS_alpha = 1e-4;
    SolverParams.LS_beta = 0.999;
    SolverParams.Minstepsize = eps;
    SolverParams.Maxstepsize = 200;
    SolverParams.Initstepsize = 1;
    SolverParams.lambdaLower = 1e-2;
    SolverParams.lambdaUpper = 1e2;
    
    SolverParams.Eps = 1e-4;
    SolverParams.Del = 1e-4;

    
    
    
    %% set general parameters for C++
    SolverParams.IsCheckParams = 0;
    SolverParams.IsCheckGradHess = 0;
    
    
    SolverParams.Max_Iteration = 200;
    SolverParams.LineSearch_LS = 4; 
    SolverParams.InitSteptype = 0;  %one
    SolverParams.isconvex = 1;
    SolverParams.HasHHR = 1;
    
    
    
    %% =========================== Compute the true solution ===========================
    TOL = 1e-7;
    Xtrue = ComputeXtrueLinfty(As, x1, TOL);
    
    
    
    %% =========================== LRBFGS-BB2 ===========================
    if (ismember('LRBFGS-BB2', TestMethods))
        for m = 1 : length(MemorySize)
            SolverParams.BBratio = 1;
            SolverParams.Num_pre_BB = 0;
            SolverParams.LengthSY = MemorySize{m};
            SolverParams.method = 'LRBFGSLPSub';            
            inFileName = ['LRBFGSLPSubBB2m' int2str(MemorySize{m})];
            TestPipeLinfty(As, x1, Xtrue, SolverParams, stabletime, SavePath, inFileName, times)
        end
    end
    
    
        
    %% =========================== NS-RBFGS ===========================
    if (ismember('RBFGS', TestMethods))
        SolverParams.method = 'RBFGSLPSub';
        inFileName = 'RBFGSLPSub';
        TestPipeLinfty(As, x1, Xtrue, SolverParams, stabletime, SavePath, inFileName, times)
    end

    
    
  %% =========================== RGS ===========================
     if (ismember('RGS', TestMethods))
            SolverParams.NumExtraGF = d + 1;
            SolverParams.method = 'RGS';
            inFileName = 'RGS';
            TestPipeLinfty(As, x1, Xtrue, SolverParams, stabletime, SavePath, inFileName, times)
    end
    
    
    %% =========================== RBFGS ===========================
    if (ismember('RBFGS-WolfeLP', TestMethods))
        SolverParams.method = 'RBFGS';
        inFileName = 'RBFGS-WolfeLP';
        TestPipeLinfty(As, x1, Xtrue, SolverParams, stabletime, SavePath, inFileName, times)
    end
    
    
    
    %% =========================== LRBFGS-WolfeLP ===========================
    if (ismember('LRBFGS-WolfeLP', TestMethods))
        for m = 1 : length(MemorySize)
            SolverParams.BBratio = 1;
            SolverParams.Num_pre_BB = 0;
            SolverParams.LengthSY = MemorySize{m};
            SolverParams.method = 'LRBFGS';            
            inFileName = ['LRBFGSBB2m' int2str(MemorySize{m}) '-WolfeLP'];
            TestPipeLinfty(As, x1, Xtrue, SolverParams, stabletime, SavePath, inFileName, times)
        end
    end
    
    
    
    %% =========================== AN ===========================
    if (ismember('AN', TestMethods))
        inFileName = 'AN';
        Max_Iteration = 5000;
        TestPipeLinftyAN(As, x1, Xtrue, Max_Iteration, stabletime, SavePath, inFileName, times)
    end

    
    
end




%% =============================================================================
% Average Data
% =============================================================================


if (ismember('RBFGS-WolfeLP', TestMethods))
    inFileName = 'RBFGS-WolfeLP';
    FileName = ['Ave' inFileName];
    for times = 1 : Ndata
        inData{times} = importdata([SavePath inFileName int2str(times) '.mat']);
        delete([SavePath inFileName num2str(times) '.mat']);
    end
     AveDataLinfty(inData, dt, SavePath, FileName);
end


if (ismember('LRBFGS-WolfeLP', TestMethods))
    for m = 1 : length(MemorySize)
        inFileName = ['LRBFGSBB2m' int2str(MemorySize{m}) '-WolfeLP'];
        FileName = ['Ave' inFileName];
        
        for times = 1 : Ndata
            inData{times} = importdata([SavePath inFileName int2str(times) '.mat']);
            delete([SavePath inFileName num2str(times) '.mat']);
        end
        AveDataLinfty(inData, dt, SavePath, FileName);
     end
end



if (ismember('RBFGS', TestMethods))
    inFileName = 'RBFGSLPSub';
    FileName = ['Ave' inFileName];
    for times = 1 : Ndata
        inData{times} = importdata([SavePath inFileName int2str(times) '.mat']);
        delete([SavePath inFileName num2str(times) '.mat']);
    end
     AveDataLinfty(inData, dt, SavePath, FileName);
end





if (ismember('LRBFGS-BB2', TestMethods))
    for m = 1 : length(MemorySize)
        inFileName = ['LRBFGSLPSubBB2m' int2str(MemorySize{m})];
        FileName = ['Ave' inFileName];
        
        for times = 1 : Ndata
            inData{times} = importdata([SavePath inFileName int2str(times) '.mat']);
            delete([SavePath inFileName num2str(times) '.mat']);
        end
        AveDataLinfty(inData, dt, SavePath, FileName);
     end
end

     




if (ismember('RGS', TestMethods))
        inFileName = 'RGS';
        FileName = ['Ave' inFileName];
        for times = 1 : Ndata
            inData{times} = importdata([SavePath inFileName int2str(times) '.mat']);
            delete([SavePath inFileName num2str(times) '.mat']);
        end
        AveDataLinfty(inData, dt, SavePath, FileName);
end




if (ismember('AN', TestMethods))
    inFileName = 'AN';
    FileName = ['Ave' inFileName];
    for times = 1 : Ndata
        inData{times} = importdata([SavePath inFileName int2str(times) '.mat']);
        %delete([SavePath inFileName num2str(times) '.mat']);
    end
    AveDataAN(inData, dt, SavePath, FileName);
end

end







