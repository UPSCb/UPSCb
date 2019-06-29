function curvature = calccurvature(X,Y,hook)
%% calccurvature - function calculates curvature base on X,Y positions
% INPUTS:
% X - x position (cm)
% Y - y position (cm)
% hook - hook distance in cm
% OUTPUT:
% curvature - sum of incremental deviations from a straight line

% Check distance from base of tree, up to hook distance. Save k for that
% distance
cSum = 0;
for k = 1:length(X)-1
  cSum = cSum +pdist([X(k:k+1),Y(k:k+1)],'euclidean');
  if(cSum>=hook), break, end
end

 xis = X(2:end) - X(1:end-1); % distance between each x (x1:x2)
 yis = Y(2:end) - Y(1:end-1); % distance between each y

%  #http://www.euclideanspace.com/maths/algebra/vectors/angleBetween/index.htm
% Calculate Angle between every pair of pixels:
for i = 1:length(xis)-1
  angle(i) = rad2deg(atan2(yis(i+1), xis(i+1))- atan2(yis(i), xis(i)));
end
% Curvature = sum of incremental deviations from a straight line:
curvature =  sum(angle);
