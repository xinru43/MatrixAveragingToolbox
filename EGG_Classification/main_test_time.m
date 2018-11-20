function [results_class, results_time] = main_test_time


%cd ToXinru;
%ImportDIR;

classes = [1, 2, 3, 4];

means = {'euclid', 'logeuclid', 'logeuclid', 'logeuclid', 'riemann', 'inductive', 'alpha', 'symmalpha'};
metric = {'euclid', 'logeuclid', 'vonbregleft', 'vonbregright', 'riemann', 'riemann', 'alpha', 'symmalpha'};

means = {'euclid', 'logeuclid', 'riemann', 'inductive', 'alpha', 'alpha', 'symmalpha', 'symmalpha'};
metric = {'euclid', 'logeuclid', 'riemann', 'riemann', 'alpha', 'alpha', 'symmalpha', 'symmalpha'};

alpha = [-1, -1, -1, -1, 0.6, 0.85, -0.9, 0.9];





N_means = length(means);


results_class = zeros(12, length(means));
results_time = zeros(12, length(means));
session = [2,2,2,2,2,2,3,2,2,4,2,5];
subjects = {'01','02','03','04','05','06','07','08','09','10','11','12'};


for sub = 1 : length(subjects)
    fprintf('--------------------------------------------------Current subject: %d\n',sub);
    for cur_session = 1 : session(sub)
        fprintf('-----session: %d\n',cur_session);
        
        % extract the batches
        s = ['EGG_Classification/data/subject' int2str(sub) '-session' int2str(cur_session) '.mat'];
        o = importdata(s);
        
        cov_train_old = o.cov_train;
        cov_test_old = o.cov_test;
        y_train = o.y_train + 1;
        y_test = o.y_test + 1;
        
        % add one empty matrices
        cov_train = cov_train_old;
        cov_test = cov_test_old;
        
        
        % -------------- loop over the algorithms computing the mean
        for i_mean = 1 : N_means
            
            % ------------------------------------------- Start training part
            fprintf(['-------Mean:' means{i_mean} '\n']);
            
            cov_centers = zeros(24, 24, length(classes));
            
            % time start: compute means of each clusters
            for i_class = 1 : length(classes)
                indx = (y_train == i_class);
                
                
                % compute the center first
                [cov_centers(:, :, i_class), ~] = mean_covariance(cov_train(:, :, indx), means{i_mean}, alpha(i_mean));
                
                
                % to get stable time
                total_time = 0;
                iter_time = 0;
                while(total_time < eps) % make sure the total running time up to 1 minute 
                    [~, each_total_time] = mean_covariance(cov_train(:, :, indx), means{i_mean}, alpha(i_mean));
                    total_time = total_time + each_total_time;
                    iter_time = iter_time + 1;
                end
                comtime(i_class) = total_time/iter_time;
                
            end
            % time end
            
            ratio = 1.0/session(sub);
            results_time(sub, i_mean) = results_time(sub, i_mean) + ratio * sum(comtime);
            % ------------------------------------------------- End training part
            
            
            
            % ------------------------------start classification on the validation set: test one by one
            for i_test = 1 : size(cov_test, 3)
                
                for i_center = 1 : size(cov_centers, 3) % loop over cluster centers
                    dist(i_center) = compute_distance(cov_test(:, :, i_test), cov_centers(:, :, i_center), metric{i_mean}, alpha(i_mean));
                end
                [~, y_pred(i_test)] = min(dist);
            end
            accuracy(cur_session, i_mean) = sum(y_pred == y_test)/length(y_test);
            % ------------------------------end classification
            
            results_class(sub, i_mean) = results_class(sub, i_mean) + ratio * accuracy(cur_session, i_mean);
        end
        % -------------- end of loop over the algorithms computing the mean
    
    end    
    
end

save('/panfs/storage.local/home-4/xy12d/ToXinru/EGG_Classification/','results_time', 'results_class', 'means', 'metric', 'alpha');

end