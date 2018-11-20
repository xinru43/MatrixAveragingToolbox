function As = L1Data(k, n, type, randomseed)

rng(randomseed);


% center = identity, not equal distance
if (type == 1)
    X = eye(n);
    L = chol(X, 'lower');
    t = 1e-1 - rand(k, 1) * 1e-2;
    for i = 1 : k
        eta{i} = randn(n);
        eta{i} = (eta{i} + eta{i}')/2;
        eta{i} = eta{i}/Metric(X, eta{i});
        As{i} = L * expm(t(k) * inv(L) * eta{i} * inv(L')) * L';
    end
end




% well conditioned data + ill conditioned outliers
if (type == 2)
    k_outlier = floor(k * 0.05);
    
    % eigenvalue mean of data points
    m = floor(n/3);
    f = 0; % for n = 3, f = 0
    mu = [rand(1, m)+1,(rand(1, n-m)+1)*10^(-f)];
    scale = 1e-3;
    
    k_original = k - k_outlier;
    for i = 1 : k_original
        % generate Gaussian vector with 0.2 standard deviation
        D = normrnd(mu, 0.2, 1, n);
        [O,~] = qr(randn(n));
        As{i} = O * diag(D) * O';
        R = chol(As{i}, 'lower');
        
        % generate white noise
        Noise = wgn(n,n,0);
        R = R + scale * Noise/norm(Noise, 2);
        As{i} = R * R';
    end
    
    % generate outliers
    m = floor(n/3);
    m = 99;
    f = 2;
    outmu = [(rand(1, n-m)+1)*10^(f), rand(1, m)+1];
    
    for i = (k_original + 1) : k
        D = normrnd(outmu, 0.2, 1, n);
        [O,~] = qr(randn(n));
        As{i} = O * diag(D) * O';
        R = chol(As{i}, 'lower');
        Noise = wgn(n,n,0);
        R = R + scale * Noise/norm(Noise, 2);
        As{i} = R * R';
    end
    
end




% ill conditioned data + well conditioned outliers
if (type == 3)
    k_outlier = floor(k * 0.05);
    
    % eigenvalue mean of data points
    if(n <= 3)
        m = 1;
    else
        m = 1;  
    end
    f = 2;
    mu = [rand(1, m)+1,(rand(1, n-m)+1)*10^(f)];
    scale = 1e-3;
    
    k_original = k - k_outlier;
    for i = 1 : k_original
        % generate Gaussian vector with 0.2 standard deviation
        D = normrnd(mu, 0.2, 1, n);
        [O,~] = qr(randn(n));
        As{i} = O * diag(D) * O';
        R = chol(As{i}, 'lower');
        
        % generate white noise
        Noise = wgn(n,n,0);
        R = R + scale * Noise/norm(Noise, 2);
        As{i} = R * R';
    end
    
    % generate outliers
    m = floor(n/3);
    f = 0;
    outmu = [(rand(1, n-m)+1)*10^(f), rand(1, m)+1];
    
    for i = (k_original + 1) : k
        D = normrnd(outmu, 0.2, 1, n);
        [O,~] = qr(randn(n));
        As{i} = O * diag(D) * O';
        R = chol(As{i}, 'lower');
        Noise = wgn(n,n,0);
        R = R + scale * Noise/norm(Noise, 2);
        As{i} = R * R';
    end
    
end




if(type == 7)
    clu = k/4;
    % cluster 1
    As{1} = eye(n);
    f = 1;
    p = 2;
    As{clu + 1} = diag([rand(1, p) + 1, (rand(1, n-p) + 1) * 10^(-f)]);
    
    p = 2;
    As{2 * clu + 1} = diag([rand(1, p) + 1, (rand(1, n-p) + 1) * 10^(f)]);
    
    f = -1;
    p = 1;
    As{3 * clu + 1} = diag([(rand(1, n-p) + 1) * 10^(f), rand(1, p) + 1]);

    for i_cluster = 1 : 4
        L = chol(As{(i_cluster - 1) * clu + 1}, 'lower');
        for j = 2 : clu
            t = 1e-2 * rand(1);
            eta = randn(n);
            eta = (eta + eta')/2;
            eta = eta/Metric(As{(i_cluster - 1) * clu + 1}, eta);
            As{(i_cluster - 1) * clu + j} = L * expm(t * inv(L) * eta * inv(L')) * L';
        end
    end
end


end



function output = Metric(X, eta)
tmp = eta * inv(X) * eta * inv(X);
output = sqrt(trace(tmp));
end