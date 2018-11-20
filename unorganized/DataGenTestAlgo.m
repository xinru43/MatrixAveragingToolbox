function As = DataGenTestAlgo(n, k, Type, TestSeed)


rng(TestSeed);
seed = rand(1, k) * 1000;
As = cell(k, 1);


% random data points
if(Type == 0)
    %f = randi([-2, 2], 1);
    f = 0;
    for i = 1 : k
        rng(seed(i));
        [O,~] = qr(randn(n));
        m = floor(n/3);
        D = [rand(1, m)+1,(rand(1, n-m)+1)*10^(f)];
        As{i} = O * diag(D) * O';
    end
end



if(Type == 1)
    for i = 1 : k
        rng(seed(i));
        l(n) = mod(i,2) + 1;
        %l(n) = i*10;
        c = (l(n))^(1.0/(n-1));
        for j = 1 : n
            l(j) = c^(j-1);
        end
        O = qr(randn(n));
        As{i} = O * diag(l) * O';
    end
end





if(Type == 2)
    L = 1;
    U = 1e+7;
    for i = 1 : k
        rng(seed(i));
        m = floor(n/2);
        s(1:m) = unifrnd(0, 0.2, [1, m]);
        s(m+1 : n) = unifrnd(0.8, 1, [1, n-m]);
        
        for j = 1 : n
            l(j) = L + (U - L)*s(n-j+1);
            
            % plot the eigenvalues
%             if (i == 1)
%                 figure(100);
%                 refline(0, l(j))
%                 hold on;
%             end
        end
        O = qr(randn(n));
        As{i} = O * diag(l) * O';
    end
end



end