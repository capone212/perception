function [C, R] = LinearPnP(X, x, K)
%% LinearPnP
% Getting pose from 2D-3D correspondences
% Inputs:
%     X - size (N x 3) matrix of 3D points
%     x - size (N x 2) matrix of 2D points whose rows correspond with X
%     K - size (3 x 3) camera calibration (intrinsics) matrix
% Outputs:
%     C - size (3 x 1) pose transation
%     R - size (3 x 1) pose rotation
%
% IMPORTANT NOTE: While theoretically you can use the x directly when solving
% for the P = [R t] matrix then use the K matrix to correct the error, this is
% more numeically unstable, and thus it is better to calibrate the x values
% before the computation of P then extract R and t directly
AA = [];
invK = inv(K);

for i = 1:size(X, 1)
    xk = invK *  [ x(i, :) 1]';
    u = xk(1);
    v = xk(2);
    Xt = [X(i,:) 1];
    A = [zeros(1,4) -Xt         v*Xt;
         Xt         zeros(1,4)  -u*Xt;
         -v*Xt      u*Xt        zeros(1,4)];
    AA = [AA; A];
end
% Solve system of linear equations using SVD
[~, ~, V] = svd(AA);
P = reshape(V(:, end), 3, 4)
[U, D, V] = svd(P(:, 1:3));
R = U*V';
C = P(:,4) / D(1,1);

if (det(U*V') < 0)
    R = -R;
    C = -C;
end

det(R)


