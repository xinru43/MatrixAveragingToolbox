function TestPipeL1EW(As, x1, Xtrue, Max_Iteration, stepsize, TOL, runingtime, SavePath, inFileName, times)

        k = length(As);
        n = size(x1, 2);


        % make data type for Bini C++
        AsEW = [];
        for i = 1 : k
            Ls{i} = chol(As{i}, 'lower');
            AsEW = [AsEW reshape(Ls{i}, 1, n * n)];
        end 
        x1EW = reshape(x1, 1, n * n);
        XtrueEW = reshape(Xtrue, 1, n * n);
    
        Debug = 2;
        [T0, ComTime0, iter, ~] = TestSPDMeanL1EW(AsEW, x1EW, XtrueEW, n, k, Max_Iteration, stepsize, TOL, Debug);
        
        T = zeros(iter+1, 1);
        
        ComTime(1) = 0;
        st = 1;
        while sum(ComTime) < runingtime
            [TT, ComTime(st)] = TestSPDMeanL1EW(AsEW, x1EW, XtrueEW, n, k, Max_Iteration, stepsize, TOL, Debug);
            T = T + TT(1:iter+1);
            st = st + 1;
        end        
        
        
        T = T/(st - 1);
        AveComTime = mean(ComTime);
        
        
        % compute distance
        Debug = 3;
        [~, ~, iter, X, diss, funs, grads, err_g] = TestSPDMeanL1EW(AsEW, x1EW, XtrueEW, n, k, Max_Iteration, stepsize, TOL, Debug);
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




