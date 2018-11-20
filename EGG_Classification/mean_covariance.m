function [mean_cov, comtime] = mean_covariance(cov_train, mean_type, alpha)


k = size(cov_train, 3);
n = size(cov_train(:, :, 1), 1);
TOL = 1e-3;


if(strcmp(mean_type, 'euclid'))
    [mean_cov, comtime] = mean_closed(cov_train, n, k, mean_type);
    mean_cov = reshape(mean_cov, n, n);
end


if(strcmp(mean_type, 'logeuclid'))
    [mean_cov, comtime] = mean_closed(cov_train, n, k, mean_type);
    mean_cov = reshape(mean_cov, n, n);
end


if(strcmp(mean_type, 'symmldbreg'))
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




end
