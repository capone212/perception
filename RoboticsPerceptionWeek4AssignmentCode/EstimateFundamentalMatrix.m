function F = EstimateFundamentalMatrix(x1, x2)
%% EstimateFundamentalMatrix
% Estimate the fundamental matrix from two image point correspondences 
% Inputs:
%     x1 - size (N x 2) matrix of points in image 1
%     x2 - size (N x 2) matrix of points in image 2, each row corresponding
%       to x1
% Output:
%    F - size (3 x 3) fundamental matrix with rank 2

A = zeros(8, 9);

for i = 1:8
    u1 = x1(i, 1);
    v1 = x1(i, 2);
    u2 = x2(i, 1);
    v2 = x2(i, 2);
    A(i, :) = [u1*u2 u1*v2 u1 v1*u2 v1*v2 v1 u2 v2 1];
end

[U D V] = svd(A);
F = reshape(V(:, 9), 3, 3)';
[U D V] = svd(F);
D(3,3) = 0;
F = U * D * V';
F = F ./norm(F);
assert(rank(F) == 2);
assert(abs(norm(F) - 1.0) < 0.0001);