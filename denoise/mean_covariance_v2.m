function [mean_cov, comtime] = mean_covariance_v2(cov_train, mean_type, alpha)

k = size(cov_train, 3);
n = size(cov_train(:, :, 1), 1);
TOL = 1e-6;


if(k == 1)
    mean_cov = cov_train(:, :, 1);
    comtime = 0;
else
    if(strcmp(mean_type, 'euclid'))
        [mean_cov, comtime] = mean_closed(cov_train, n, k, mean_type);
        mean_cov = reshape(mean_cov, n, n);
    end
    
    if(strcmp(mean_type, 'logeuclid'))
        [mean_cov, comtime] = mean_closed(cov_train, n, k, mean_type);
        mean_cov = reshape(mean_cov, n, n);
    end
    
    
    
    if(strcmp(mean_type, 'inductive'))
        [mean_cov, comtime] = mean_closed(cov_train, n, k, mean_type);
        mean_cov = reshape(mean_cov, n, n);
    end
    
    
    
    if(strcmp(mean_type, 'riemann'))
        [mean_cov, comtime] = mean_riemann(cov_train, TOL);
    end
    
    
    
    if(strcmp(mean_type, 'alpha'))
        [mean_cov, comtime] = mean_alpha(cov_train, alpha, TOL);
    end
    
    
    
    if(strcmp(mean_type, 'symmalpha'))
        [mean_cov, comtime] = mean_symmalpha(cov_train, alpha, TOL);
    end
    
    
    if(strcmp(mean_type, 'symmldbreg'))
        [mean_cov, comtime] = mean_closed(cov_train, n, k, mean_type);
        mean_cov = reshape(mean_cov, n, n);
    end
    
    
    % meidan
    
    if(strcmp(mean_type, 'medianriemann'))
        [mean_cov, comtime] = median_riemann(cov_train, TOL);
    end
    
    
    if(strcmp(mean_type, 'medianalpha'))
        [mean_cov, comtime] = median_alpha(cov_train, alpha, TOL);
    end
    
    
    if(strcmp(mean_type, 'mediansymmalpha'))
        [mean_cov, comtime] = median_symmalpha(cov_train, alpha, TOL);
    end
    

    % circumcenter   
    if(strcmp(mean_type, 'minimaxcenterriemann'))
        [mean_cov, comtime] = minimaxcenter_riemann(cov_train, TOL);
    end
                            
    if(strcmp(mean_type, 'minimaxcenteralpha'))
        [mean_cov, comtime] = minimaxcenter_alpha(cov_train, alpha, TOL);
    end
    
    if(strcmp(mean_type, 'minimaxcentersymmalpha'))
        [mean_cov, comtime] = minimaxcenter_symmalpha(cov_train, alpha, TOL);
    end
    
    
    
end

% if(strcmp(mean_type, 'harmonic'))
%     mean_cov = inv(cov_train(:, : , 1));
%     for i = 2 : k
%         mean_cov = mean_cov + inv(cov_train(:, : , i));
%     end
%     mean_cov = pinv(mean_cov/k);
%     comtime = 1;
% end


%     Ar = cov_train(:, : , 1);
%     for i = 2 : k
%         Ar = Ar + cov_train(:, : , i);
%     end
%     Ar = Ar/k;
%
%     Ha = inv(cov_train(:, : , 1));
%     for i = 2 : k
%         Ha = Ha + inv(cov_train(:, : , i));
%     end
%     Ha = inv(Ha/k);
%
%     Arh = Ar^(0.5);
%     Arinvh = Ar^(-0.5);
%     temp = Arinvh * Ha * Arinvh;
%     mean_cov = Arh * temp^(0.5) * Arh;
%     here = norm(mean_cov - mean_cov2)
end
