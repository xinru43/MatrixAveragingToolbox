function TestPipeL2RSD(As, x1, Xtrue, Max_Iteration, Tolerance, Retraction, StepType, runingtime, SavePath, inFileName, times)

        k = length(As);
        n = size(x1, 2);


        % make data type for Bini C++
        As_ = [];
        for i = 1 : k
            Ls{i} = chol(As{i}, 'lower');
            As_ = [As_ reshape(Ls{i}, 1, n * n)];
        end 
        x1_ = reshape(x1, 1, n * n);
        Xtrue_ = reshape(Xtrue, 1, n * n);
    
        Debug = 2;
        [T0, ComTime0, iter, ~] = TestSPDMeanRSD(As_, x1_, Xtrue_, n, k, Max_Iteration, Tolerance, Retraction, StepType, Debug);
        
        T = zeros(iter, 1);
        
        ComTime(1) = 0;
        st = 1;
        while sum(ComTime) < runingtime
            [TT, ComTime(st)] = TestSPDMeanRSD(As_, x1_, Xtrue_, n, k, Max_Iteration, Tolerance, Retraction, StepType, Debug);
            T = T + TT(1:iter);
            st = st + 1;
        end        
        
        
        T = T/(st - 1);
        AveComTime = mean(ComTime);
        
        % compute distance
        Debug = 3;
        [~, ~, iter, X, diss, funs, grads, err_g] = TestSPDMeanRSD(As_, x1_, Xtrue_, n, k, Max_Iteration, Tolerance, Retraction, StepType, Debug);
        dis = diss(1:iter);
        G = grads(1:iter);
        F = funs(1:iter);
        
        T = T';
        dis = dis';
        F = F';
        G = G';

        save([SavePath inFileName int2str(times) '.mat'],'ComTime', 'AveComTime',...
        'G','T', 'F','dis','iter', 'T0', 'ComTime0');

end




