function output = compute_distance(x, y, dist_type, alpha)
% x = input
% y = cluster center


if(strcmp(dist_type, 'euclid'))
    temp = x - y;
    output = trace(temp' * temp);
end


if(strcmp(dist_type, 'logeuclid'))
    temp = logm(x) - logm(y);
    output = trace(temp' * temp);
end



if(strcmp(dist_type, 'riemann')) 
    L = chol(x, 'lower');
    temp = (L\y)/L';
    temp = logm(temp);
    output = trace(temp' * temp);
end



if(strcmp(dist_type, 'inductive'))
    L = chol(x, 'lower');
    temp = (L\b)/L';
    temp = logm(temp);
    output = trace(temp' * temp);
end


if(strcmp(dist_type, 'alpha'))
    a = (1.0 - alpha)/2.0;
    b = (1.0 + alpha)/2.0;
    output = log(det(a * x + b * y)) - a * log(det(x)) - b * log(det(y));
end



if(strcmp(dist_type, 'symmalpha'))
    a = (1.0 - alpha)/2.0;
    b = (1.0 + alpha)/2.0;
    output = log(det(a * x + b * y)) + log(det(b * x + a * y)) - log(det(x)) - log(det(y));
end



if(strcmp(dist_type, 'symmldbreg'))
    output = trace(y\x + x\y - 2 * eye(size(x)));
end



% delta(A, mean)
if(strcmp(dist_type, 'vonbregleft'))
    temp = x * logm(x) - x * logm(y) - x + y;
    output = trace(temp);
end


% delta(mean, A)
if(strcmp(dist_type, 'vonbregright'))
    temp = y * logm(y) - y * logm(x) - y + x;
    output = trace(temp);
end


% symmetric von neumann bregman
if(strcmp(dist_type, 'vonbregsymm'))
    temp = y * logm(y) - y * logm(x) - y + x;
    output = trace(temp);
    
    temp = x * logm(x) - x * logm(y) - x + y;
    output = output + trace(temp);
    output = output/2.0;
end


end