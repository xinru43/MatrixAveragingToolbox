function main_denoise

n = 128;
W = 3;


T = importdata('/gpfs/home/xy12d/ToXinru/average/duorou3_T.mat');
T_hat = importdata('/gpfs/home/xy12d/ToXinru/average/duorou3_T_hat_pr002.mat');

%% denoise

for i = 1 : 10
    
fprintf('=========== iter: %d\n', i); 
mean_type = 'euclid'
[dT_euclid{i}, time_euclid{i}] = denoise(T_hat{i}, W, mean_type, 0);

mean_type = 'logeuclid'
[dT_logeuclid{i}, time_logeuclid{i}] = denoise(T_hat{i}, W, mean_type, 0);

mean_type = 'symmldbreg'
[dT_symmldbreg{i}, time_symmldbreg{i}] = denoise(T_hat{i}, W, mean_type, 0);


mean_type = 'riemann'
[dT_riemann{i}, time_riemann{i}] = denoise(T_hat{i}, W, mean_type, 0);

mean_type = 'medianriemann'
[dT_medianriemann{i}, time_medianriemann{i}] = denoise(T_hat{i}, W, mean_type, 0);


mean_type = 'minimaxcenterriemann'
[dT_minimaxcenterriemann{i}, time_minimaxcenterriemann{i}] = denoise(T_hat{i}, W, mean_type, 0);
 

alpha = 0
mean_type = 'alpha'
[dT_alpha0{i}, time_alpha0{i}] = denoise(T_hat{i}, W, mean_type, alpha);

mean_type = 'medianalpha'
[dT_medianalpha0{i}, time_medianalpha0{i}] = denoise(T_hat{i}, W, mean_type, alpha);

mean_type = 'minimaxcenteralpha'
[dT_minimaxcenteralpha0{i}, time_minimaxcenteralpha0{i}] = denoise(T_hat{i}, W, mean_type, alpha);



alpha = 0.1
mean_type = 'alpha'
[dT_alpha1{i}, time_alpha1{i}] = denoise(T_hat{i}, W, mean_type, alpha);

mean_type = 'medianalpha'
[dT_medianalpha1{i}, time_medianalpha1{i}] = denoise(T_hat{i}, W, mean_type, alpha);

mean_type = 'minimaxcenteralpha'
[dT_minimaxcenteralpha1{i}, time_minimaxcenteralpha1{i}] = denoise(T_hat{i}, W, mean_type, alpha);




alpha = 0.2
mean_type = 'alpha'
[dT_alpha2{i}, time_alpha2{i}] = denoise(T_hat{i}, W, mean_type, alpha);

mean_type = 'medianalpha'
[dT_medianalpha2{i}, time_medianalpha2{i}] = denoise(T_hat{i}, W, mean_type, alpha);

mean_type = 'minimaxcenteralpha'
[dT_minimaxcenteralpha2{i}, time_minimaxcenteralpha2{i}] = denoise(T_hat{i}, W, mean_type, alpha);



alpha = 0.3
mean_type = 'alpha'
[dT_alpha3{i}, time_alpha3{i}] = denoise(T_hat{i}, W, mean_type, alpha);

mean_type = 'medianalpha'
[dT_medianalpha3{i}, time_medianalpha3{i}] = denoise(T_hat{i}, W, mean_type, alpha);

mean_type = 'minimaxcenteralpha'
[dT_minimaxcenteralpha3{i}, time_minimaxcenteralpha3{i}] = denoise(T_hat{i}, W, mean_type, alpha);




alpha = 0.4
mean_type = 'alpha'
[dT_alpha4{i}, time_alpha4{i}] = denoise(T_hat{i}, W, mean_type, alpha);

mean_type = 'medianalpha'
[dT_medianalpha4{i}, time_medianalpha4{i}] = denoise(T_hat{i}, W, mean_type, alpha);

mean_type = 'minimaxcenteralpha'
[dT_minimaxcenteralpha4{i}, time_minimaxcenteralpha4{i}] = denoise(T_hat{i}, W, mean_type, alpha);



alpha = 0.5
mean_type = 'alpha'
[dT_alpha5{i}, time_alpha5{i}] = denoise(T_hat{i}, W, mean_type, alpha);

mean_type = 'medianalpha'
[dT_medianalpha5{i}, time_medianalpha5{i}] = denoise(T_hat{i}, W, mean_type, alpha);

mean_type = 'minimaxcenteralpha'
[dT_minimaxcenteralpha5{i}, time_minimaxcenteralpha5{i}] = denoise(T_hat{i}, W, mean_type, alpha);



alpha = 0.6
mean_type = 'alpha'
[dT_alpha6{i}, time_alpha6{i}] = denoise(T_hat{i}, W, mean_type, alpha);

mean_type = 'medianalpha'
[dT_medianalpha6{i}, time_medianalpha6{i}] = denoise(T_hat{i}, W, mean_type, alpha);

mean_type = 'minimaxcenteralpha'
[dT_minimaxcenteralpha6{i}, time_minimaxcenteralpha6{i}] = denoise(T_hat{i}, W, mean_type, alpha);



alpha = 0.7
mean_type = 'alpha'
[dT_alpha7{i}, time_alpha7{i}] = denoise(T_hat{i}, W, mean_type, alpha);

mean_type = 'medianalpha'
[dT_medianalpha7{i}, time_medianalpha7{i}] = denoise(T_hat{i}, W, mean_type, alpha);

mean_type = 'minimaxcenteralpha'
[dT_minimaxcenteralpha7{i}, time_minimaxcenteralpha7{i}] = denoise(T_hat{i}, W, mean_type, alpha);



alpha = 0.8
mean_type = 'alpha'
[dT_alpha8{i}, time_alpha8{i}] = denoise(T_hat{i}, W, mean_type, alpha);

mean_type = 'medianalpha'
[dT_medianalpha8{i}, time_medianalpha8{i}] = denoise(T_hat{i}, W, mean_type, alpha);

mean_type = 'minimaxcenteralpha'
[dT_minimaxcenteralpha8{i}, time_minimaxcenteralpha8{i}] = denoise(T_hat{i}, W, mean_type, alpha);



alpha = 0.9
mean_type = 'alpha'
[dT_alpha9{i}, time_alpha9{i}] = denoise(T_hat{i}, W, mean_type, alpha);

mean_type = 'medianalpha'
[dT_medianalpha9{i}, time_medianalpha9{i}] = denoise(T_hat{i}, W, mean_type, alpha);

mean_type = 'minimaxcenteralpha'
[dT_minimaxcenteralpha9{i}, time_minimaxcenteralpha9{i}] = denoise(T_hat{i}, W, mean_type, alpha);


end

 save('/gpfs/home/xy12d/ToXinru/average/denoise_time_pr002.mat', ...
 'time_euclid', 'time_logeuclid', 'time_symmldbreg', ...
 'time_riemann', 'time_medianriemann', 'time_minimaxcenterriemann', ...
 'time_alpha0','time_medianalpha0', 'time_minimaxcenteralpha0', ...
 'time_alpha1','time_medianalpha1', 'time_minimaxcenteralpha1', ...
 'time_alpha2','time_medianalpha2', 'time_minimaxcenteralpha2', ...
 'time_alpha3','time_medianalpha3', 'time_minimaxcenteralpha3', ...
 'time_alpha4','time_medianalpha4', 'time_minimaxcenteralpha4', ...
 'time_alpha5','time_medianalpha5', 'time_minimaxcenteralpha5', ...
 'time_alpha6','time_medianalpha6', 'time_minimaxcenteralpha6', ...
 'time_alpha7','time_medianalpha7', 'time_minimaxcenteralpha7', ...
 'time_alpha8','time_medianalpha8', 'time_minimaxcenteralpha8', ...
 'time_alpha9','time_medianalpha9', 'time_minimaxcenteralpha9');
 
 
  save('/gpfs/home/xy12d/ToXinru/average/denoise_result_pr002.mat', ...
 'dT_euclid', 'dT_logeuclid', 'dT_symmldbreg', ...
 'dT_riemann', 'dT_medianriemann', 'dT_minimaxcenterriemann', ...
 'dT_alpha0','dT_medianalpha0', 'dT_minimaxcenteralpha0', ...
 'dT_alpha1','dT_medianalpha1', 'dT_minimaxcenteralpha1', ...
 'dT_alpha2','dT_medianalpha2', 'dT_minimaxcenteralpha2', ...
 'dT_alpha3','dT_medianalpha3', 'dT_minimaxcenteralpha3', ...
 'dT_alpha4','dT_medianalpha4', 'dT_minimaxcenteralpha4', ...
 'dT_alpha5','dT_medianalpha5', 'dT_minimaxcenteralpha5', ...
 'dT_alpha6','dT_medianalpha6', 'dT_minimaxcenteralpha6', ...
 'dT_alpha7','dT_medianalpha7', 'dT_minimaxcenteralpha7', ...
 'dT_alpha8','dT_medianalpha8', 'dT_minimaxcenteralpha8', ...
 'dT_alpha9','dT_medianalpha9', 'dT_minimaxcenteralpha9');
 

end



function [dT, comtime] = denoise(T, W, mean_type, alpha)

[tm, tn, m, n] = size(T);
dT = zeros(tm, tn, m, n);

progress_flag = 0;

for i = 1 : tm
    for j = 1 : tn
        % get the neighbor index of pixel (i, j)
        neighbor = get_window_idx(i, j, tm, tn, W);
        
        % get the neighbor tensors and stored in As
        for k = 1 : length(neighbor)
            ni = neighbor{k}(1);
            nj = neighbor{k}(2);
            As(:, :, k) = reshape(T(ni, nj, :, :), m, n);
        end
        % compute the average in the neighbor
         [dT(i, j, :, :), comtime(i, j)] = mean_covariance_v2(As, mean_type, alpha);
         progress_flag = progress_flag + 1;
         
         if(progress_flag == 1638)
             fprintf('=========== 10%% \n');
         elseif(progress_flag == 1638 * 2)
             fprintf('=========== 20%% \n');
         elseif(progress_flag == 1638 * 3)
             fprintf('=========== 30%% \n');
         elseif(progress_flag == 1638 * 4)
             fprintf('=========== 40%% \n');
         elseif(progress_flag == 1638 * 5)
             fprintf('=========== 50%% \n');
         elseif(progress_flag == 1638 * 6)
             fprintf('=========== 60%% \n');
         elseif(progress_flag == 1638 * 7)
             fprintf('=========== 70%% \n');
         elseif(progress_flag == 1638 * 8)
             fprintf('=========== 80%% \n');
         elseif(progress_flag == 1638 * 9)
             fprintf('=========== 90%% \n');
         elseif(progress_flag == 1638 * 10)
             fprintf('=========== 100%% \n');
         end
         
    end
end

end



function output = check_noise(neighbor, noise_idx)
% check if this window contains noise
output = 0;
for k = 1 : length(neighbor)
    i = neighbor{k}(1);
    j = neighbor{k}(2);
    neighbor_idx = [int2str(i) '&' int2str(j)];
    tmp = ismember(neighbor_idx, noise_idx);
    if(tmp == 1)
        output = 1;
        break;
    end
end

end







function output = get_window_idx(i, j, m, n, w)

b = (w-1)/2;
i_idx = [i - b : i + b];
j_idx = [j - b : j + b];

i_idx = i_idx(i_idx >= 1);
i_idx = i_idx(i_idx <= m);

j_idx = j_idx(j_idx >= 1);
j_idx = j_idx(j_idx <= n);


k = 1;
for a = 1 : length(i_idx)
    for b = 1 : length(j_idx)
        output{k} = [i_idx(a) j_idx(b)];
        k = k + 1;
    end
end


end







% % check positive definite of tensors
% cumm = 0;
% for i = 1 : n
%     for j = 1 : n
%         a = reshape(T(i, j, :, :), 2, 2);
%         kappa(i, j) = cond(a);
%         l = eig(a);
%         tmp = sum(l > 0);
%         if (tmp < 2)
%             cumm = cumm + 1;
%         end
%     end
% end



