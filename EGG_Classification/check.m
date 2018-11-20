n = 50;
A = randn(n);
B = randn(n);
A = A * A';
B = B * B';

t = 1.0/3;
v1 = A^(0.5) * (A^(-0.5) * B * A^(-0.5))^(t) * A^(0.5);

a = chol(A, 'lower');
b = chol(B, 'lower');
v2 = a * (inv(a) * B * inv(a'))^t * a';

norm(v1 -  v2)



n = 3;
k = 2;
A = zeros(n, n, k);
for i = 1 : k
    A(:, :, i)=randn(n);
    A(:, :, i) = A(:, :, i) * A(:, :, i)';
end


n = 3;
k = 2;
t = 1.0/2;
a = chol(A, 'lower');
b = chol(B, 'lower');
Z = inv(a) * b;
C = Z * Z';
(Z * Z')^t



