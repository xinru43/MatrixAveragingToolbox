function TestPipeLDOneParamLinftyAN(As, alpha, x1, Xtrue, Max_Iteration, runingtime, SavePath, inFileName, times)


        k = length(As);
        n = size(x1, 2);
        TOL = 1e-20;


        % make data type for Bini C++
        AsAN = [];
        for i = 1 : k
            Ls{i} = chol(As{i}, 'lower');
            AsAN = [AsAN reshape(Ls{i}, 1, n * n)];
            LDAs(i) = log(det(As{i}));
        end 
        x1AN = reshape(x1, 1, n * n);
        XtrueAN = reshape(Xtrue, 1, n * n);
    
        Debug = 2;
        [T0, ComTime0, iter, ~] = TestSPDMeanLDOneParamLinftyAN(AsAN, LDAs, alpha, x1AN, XtrueAN, n, k, Max_Iteration, 1e-7, Debug);
        
        T = zeros(iter+1, 1);
        
        ComTime(1) = 0;
        st = 1;
        while sum(ComTime) < runingtime
            [TT, ComTime(st)] = TestSPDMeanLDOneParamLinftyAN(AsAN, LDAs, alpha, x1AN, XtrueAN, n, k, Max_Iteration, 1e-7, Debug);
            T = T + TT(1:iter+1);
            st = st + 1;
        end        
        
        
        T = T/(st - 1);
        AveComTime = mean(ComTime);
        
        
        % compute distance
        Debug = 3;
        [~, ~, iter, X, diss, funs, grads, err_g] = TestSPDMeanLDOneParamLinftyAN(AsAN, LDAs, alpha, x1AN, XtrueAN, n, k, Max_Iteration, 1e-7, Debug);
        dis = diss(1:iter+1);
        G = grads(1:iter+1);
        F = funs(1:iter+1);
        
        T = T';
        dis = dis';
        F = F';
        G = G';

        save([SavePath inFileName int2str(times) '.mat'],'ComTime', 'AveComTime',...
        'G','T', 'F','dis','iter', 'T0', 'ComTime0');

end




