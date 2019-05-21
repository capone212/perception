function [ H ] = est_homography(video_pts, logo_pts)
% est_homography estimates the homography to transform each of the
% video_pts into the logo_pts
% Inputs:
%     video_pts: a 4x2 matrix of corner points in the video
%     logo_pts: a 4x2 matrix of logo points that correspond to video_pts
% Outputs:
%     H: a 3x3 homography matrix such that logo_pts ~ H*video_pts
% Written for the University of Pennsylvania's Robotics:Perception course

% YOUR CODE HERE
A = [];
for i = 1:4
  x = video_pts(i, :);
  xs = logo_pts(i, :);
  ax = [-x(1), -x(2), -1, 0, 0, 0, x(1)*xs(1), x(2)*xs(1), xs(1)];
  ay = [0, 0, 0, -x(1), -x(2), -1, x(1)*xs(2), x(2)*xs(2), xs(2)];
  A = [A; ax; ay];
end
[U, S, V] = svd(A);

h = V(:,9)';

H = [h(1:3); h(4:6); h(7:9)];

end

