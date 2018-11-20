function As = GenerateDataToy(n, k, ko)
% Generate two classes
seed = 1000;
rng(seed);

tol = 1.0e+1;


% the eigenvalue distribution of data points
% 
% for i = 1 : k
%     
%     D = diag([1 + 5 * rand(n,1)]);
%     [O,~] = qr(randn(n));
%     As{i} = O * D * O';
%     
% % %     [O,~] = qr(randn(n));
% % %     D = diag([rand(1,n-m)+1,(rand(1,m)+1)*10^(-f)]);
% % %     As{i} = O * D * O';
% %     C = diag([rand(1,n-m)+1,(rand(1,m)+1)*10^(f)]);
% %     % noise 
% %     D = diag(rand(1,n)+1);
% %     [O,~] = qr(randn(n));
% %     noise = O * D * O';
% %     
% % %     [O,~] = qr(randn(n));
% %     As{i} = O * C * O' + tol * noise/norm(noise, 2);
% end




L = 1;
U = 2;
m = floor(n/2);
for i = 1 : k
    s(1:m) = unifrnd(0, 0.4, [1, m]);
    s(m+1 : n) = unifrnd(0.8, 1, [1, n-m]);
    
    for j = 1 : n
        l(j) = L + (U - L)*s(n-j+1);
                    % plot the eigenvalues
                    if (i == 1)
                        figure(100);
                        refline(0, l(j))
                        hold on;
                    end
    end
    O = qr(randn(n));
    As{i} = O * diag(l) * O';
end



% fprintf('Within class 1 distance: \n');
% for i = 1 : k
%     for j = i+1 : k
%         WC1(i,j) = AIRM(As{i}, As{j})
%     end
% end

for i = 1 : ko
D = diag([5*(rand(1,n - m)+1)*10^(6), (rand(1,m)+1) * 10^(6)]);
[O,~] = qr(randn(n));
As{k + i} = O * D * O';
end

% fprintf('Distance from outliers: \n');
% for i = 1 : k
%     OutlierToC1(i) = AIRM(As{k + 1}, As{i})
% end



end