function As = LinftyData(k, n, type, randomseed)

rng(randomseed);

% center = identity, equal distance, distance = 1e-2
if (type == 1)
    t = 1e-1;
    X = eye(n);
    L = chol(X, 'lower');
    for i = 1 : k
        eta{i} = randn(n);
        eta{i} = (eta{i} + eta{i}')/2;
        eta{i} = eta{i}/Metric(X, eta{i});
        As{i} = L * expm(t * inv(L) * eta{i} * inv(L')) * L';
    end
end


% center = identity, not equal distance
if (type == 2)
    X = eye(n);
    L = chol(X, 'lower');
    %t = 10 * rand(k, 1) + 1e-2;
    t = 1e-1 - rand(k, 1) * 1e-3;
    for i = 1 : k
        eta{i} = randn(n);
        eta{i} = (eta{i} + eta{i}')/2;
        eta{i} = eta{i}/Metric(X, eta{i});
        As{i} = L * expm(t(k) * inv(L) * eta{i} * inv(L')) * L';
    end
end


% center = ill conditioned, equal distance, distance = 1
if(type == 3)
    t = 1;
    p = floor(n/3);
    f = 5;
    [O, ~] = qr(randn(n));
    D = diag([rand(1, p) + 1, (rand(1, n-p) + 1) * 10^(-f)]);
    X = O * D * O';
    L = chol(X, 'lower');
    for i = 1 : k
        eta{i} = randn(n);
        eta{i} = (eta{i} + eta{i}')/2;
        eta{i} = eta{i}/Metric(X, eta{i});
        As{i} = L * expm(t * inv(L) * eta{i} * inv(L')) * L';
    end
end




% center = ill conditioned, not equal distance
if(type == 4)
    p = floor(n/3);
    f = 5;
    [O, ~] = qr(randn(n));
    D = diag([rand(1, p) + 1, (rand(1, n-p) + 1) * 10^(-f)]);
    X = O * D * O';
    L = chol(X, 'lower');
    
    t = 1 + rand(k, 1);

    for i = 1 : k
        eta{i} = randn(n);
        eta{i} = (eta{i} + eta{i}')/2;
        eta{i} = eta{i}/Metric(X, eta{i});
        As{i} = L * expm(t(i) * inv(L) * eta{i} * inv(L')) * L';
    end
end



% a few clusters
if(type == 5)
   
    t = 1e-4 * rand(k, 1);
    
    % cluster 1
    As{1} = eye(n);
    L = chol(As{1}, 'lower');
    for i = 2 : 5
        eta = randn(n);
        eta = (eta + eta')/2;
        eta = eta/Metric(As{1}, eta); 
        As{i} = L * expm(t(i) * inv(L) * eta * inv(L')) * L';
    end
    
    
    % cluster 2
    L = chol(As{1}, 'lower');
    d = randn(n);
    d = (d + d')/2;
    d = d/Metric(As{1}, d);
    scale = 1e+1;
    As{6} = L * expm(scale * inv(L) * d * inv(L')) * L';
    L = chol(As{6}, 'lower');
    for i = 7 : 10
        eta = randn(n);
        eta = (eta + eta')/2;
        eta = eta/Metric(As{6}, eta);
        As{i} = L * expm(t(i) * inv(L) * eta * inv(L')) * L';
    end
    
    
    % cluster 3
    L = chol(As{6}, 'lower');
    d = randn(n);
    d = (d + d')/2;
    d = d/Metric(As{6}, d);
    scale = 1e+1;
    As{11} = L * expm(scale * inv(L) * d * inv(L')) * L';
    L = chol(As{11}, 'lower');
    for i = 12 : 15
        eta = randn(n);
        eta = (eta + eta')/2;
        eta = eta/Metric(As{11}, eta);
        As{i} = L * expm(t(i) * inv(L) * eta * inv(L')) * L';
    end
    
    
    % cluster 4
    L = chol(As{11}, 'lower');
    d = randn(n);d = (d + d')/2;
    d = d/Metric(As{11}, d);
    scale = 5;
    As{16} = L * expm(scale * inv(L) * d * inv(L')) * L';
    L = chol(As{16}, 'lower');
    for i = 17 : 20
        eta = randn(n);
        eta = (eta + eta')/2;
        eta = eta/Metric(As{16}, eta);
        As{i} = L * expm(t(i) * inv(L) * eta * inv(L')) * L';
    end
    
end




if(type == 6)
    clu = k/4;
    % cluster 1
    As{1} = eye(n);
    f = 4;
    p = 1;
    As{clu + 1} = diag([rand(1, p) + 1, (rand(1, n-p) + 1) * 10^(-f)]);
    
    p = 2;
    As{2 * clu + 1} = diag([rand(1, p) + 1, (rand(1, n-p) + 1) * 10^(f)]);
    
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



if(type == 7)
    clu = k/4;
    % cluster 1
    As{1} = eye(n);
    f = 2;
    p = 1;
    As{clu + 1} = diag([rand(1, p) + 1, (rand(1, n-p) + 1) * 10^(-f)]);
    
    p = 2;
    As{2 * clu + 1} = diag([rand(1, p) + 1, (rand(1, n-p) + 1) * 10^(f)]);
    
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