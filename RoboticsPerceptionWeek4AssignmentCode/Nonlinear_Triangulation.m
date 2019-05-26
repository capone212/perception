function X = Nonlinear_Triangulation(K, C1, R1, C2, R2, C3, R3, x1, x2, x3, X0)
%% Nonlinear_Triangulation
% Refining the poses of the cameras to get a better estimate of the points
% 3D position
% Inputs: 
%     K - size (3 x 3) camera calibration (intrinsics) matrix
%     x
% Outputs: 
%     X - size (N x 3) matrix of refined point 3D locations 
    %XR = Reproject(C1, R1, K, X0(1, :)');
    X = [];
    for i = 1:size(X0,1)
        XS = Single_Point_Nonlinear_Triangulation(K, C1, R1, C2, R2, C3, R3, x1(i,:)', x2(i,:)', x3(i,:)', X0(i,:)');
        X = [X;XS'];
    end

end

function dX = deltaX(K, C1, R1, C2, R2, C3, R3, x1, x2, x3, X)
    XR1 = Reproject(C1, R1, K, X);
    XR2 = Reproject(C2, R2, K, X);
    XR3 = Reproject(C3, R3, K, X);
    F = FX(XR1, XR2, XR3);
    J = Jacobian(K, XR1, XR2,XR3, R1, R2, R3);
    b = [x1(1) x1(2) x2(1) x2(2) x3(1) x3(2)]';
    dX = inv(J'*J)*J'*(b-F);
end

function X = Single_Point_Nonlinear_Triangulation(K, C1, R1, C2, R2, C3, R3, x1, x2, x3, X0)
    X = X0;
    for i = 1:3
        dx = deltaX(K, C1, R1, C2, R2, C3, R3, x1, x2, x3, X);
        X = X + dx;
    end
end


function J = Jacobian(K, XR1, XR2, XR3, R1, R2, R3)
    J = [dF_dX(XR1, K,  R1)' dF_dX(XR2, K,  R2)'   dF_dX(XR3, K,  R3)']';
end

function F = FX(XR1, XR2, XR3)
    F = [scale_XR(XR1) scale_XR(XR2) scale_XR(XR3)]';
end

function XR = Reproject(C, R, K, X)
    XR = K * R * (X - C);
end

function S = scale_XR(XR)
    S = [XR(1)/XR(3) XR(2)/XR(3)];
end

function dfdx = dF_dX(XR, K,  R)
    u = XR(1);
    v = XR(2);
    w = XR(3);
    part = K*R;
    dudx = part(1,:);
    dvdx = part(2,:);
    dwdx = part(3,:);
    dfdx = [(w*dudx - u*dwdx)/w^2;
            (w*dvdx - v*dwdx)/w^2];
end

function dfdx = dF_dX2(FX, K,  R)
    u = FX(1);
    v = FX(2);
    w = FX(3);
    dudx = du_dX(K, R);
    dvdx = dv_dX(K, R);
    dwdx = dw_dx(R);
    dfdx = [(w*dudx - u*dwdx)/w^2;
            (w*dvdx - v*dwdx)/w^2];
end

function dudx = du_dX(K, R)
    [f, px, py] = f_px_py(K);
    dudx = [f*R(1,1)+px*R(3,1) f*R(1,2)+px*R(3,2) f*R(1,3)+px*R(3,3)];
end

function dvdx = dv_dX(K, R)
    [f, px, py] = f_px_py(K);
    dvdx = [f*R(2,1)+py*R(3,1) f*R(2,2)+py*R(3,2) f*R(2,3)+py*R(3,3)];
end

function dwdx = dw_dx(R)
    dwdx = R(3, :);
end

function [f, px, py] = f_px_py(K)
    f = K(1,1);
    px = K(1, 3);
    py = K(2, 3);
end
