function TestPipeFP(As, alpha, x1, Xtrue, Max_Iteration, stabletime, SavePath, inFileName, times)


        k = length(As);
        n = size(x1, 2);
        TOL = 1e-20;


        % make data type for Bini C++
        AsFP = [];
        for i = 1 : k
            AsFP = [AsFP reshape(As{i}, 1, n * n)];
        end 
        x1FP = reshape(x1, 1, n * n);
        XtrueFP = reshape(Xtrue, 1, n * n);
    
        Debug = 2;
        [T0, ComTime0, iter, ~] = TestSPDDLOneParamFP(AsFP, alpha, x1FP, XtrueFP, n, k, Max_Iteration, 1e-7, Debug);
        
        T = zeros(iter+1, 1);
        for st = 1 : stabletime
            [TT, ComTime(st), iter, ~] = TestSPDDLOneParamFP(AsFP, alpha, x1FP, XtrueFP, n, k, Max_Iteration, 1e-7, Debug);
            T = T + TT(1:iter+1);
        end
        T = T/stabletime;
        AveComTime = mean(ComTime);
        
        
        % compute distance
        Debug = 3;
        [~, ~, iter, X, diss, truestepsizes, Gs, gfgf0] = TestSPDDLOneParamFP(AsFP, alpha, x1FP, XtrueFP, n, k, Max_Iteration, 1e-7, Debug);
        dis = diss(1:iter+1);
        truestepsize = truestepsizes(1 : iter+1);
        G = Gs(1:iter+1);
        T = T';
        dis = dis';
        G = G';

        save([SavePath inFileName int2str(times) '.mat'],'ComTime', 'AveComTime',...
        'G','T','dis','iter', 'truestepsize', 'T0', 'ComTime0', 'gfgf0');

end




