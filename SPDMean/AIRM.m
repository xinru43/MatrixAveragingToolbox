function output = AIRM(X, Y)  % compute distance


Lx = chol(X, 'lower');
D = eig(Lx\Y/Lx');
output = sum(log(D).^2);
output = sqrt(output);

% Xinvh = X^(-0.5);
% [V, D] = schur(Xinvh * Y * Xinvh);
% output = sum(log(diag(D)).^2);
% output = sqrt(output)
% here = output - output2
end

