function [x1_curr, x2_curr, T] = Ransac(x1, x2)

numIter = 2000;
best = 0;
thresh = 49;

% Run RANSAC for required number of iterations
for i=1:numIter
  % Select a subset of points (4)
  samples = randi([1 size(x1, 1)], 4, 1);
  x1_points = x1(samples, :);
  x2_points = x2(samples, :);

  % Using the selected points, estimate homography
  p = 1;
  for c = 1:size(x1_points,1)
    x_h = x1_points(c,1); 
    xprime_h = x2_points(c,1); 
    y_h = x1_points(c,2);
    yprime_h = x2_points(c,2);
    %use homography matrix structure from lecture slides
    A(p,:) = [x_h y_h 1 0 0 0 -xprime_h*x_h -xprime_h*y_h -xprime_h];
    A(p+1,:) = [0 0 0 x_h y_h 1 -yprime_h*x_h -yprime_h*y_h -yprime_h];

    p = p+2;
  end
  
  % must find A'A with smallest eigenvalue (need diagonalized matrix) to
  % get homography matrix
  [~,~,H] = svd(A);
  
  %reshape homography matrix to be 3x3
  H = reshape(H(:,end),3,3);
  %make projective transformation
  t = maketform('projective',H);
  [x1_transform(:,1) x1_transform(:,2)] = tformfwd(t,x1(:,1),x1(:,2));
  
  % Score homography by counting number of points
  % lying within the threshold
  val = sum(abs(x2 - x1_transform).^2,2);
  curr = val < thresh;
  
  % If this is the highest-scoring homography found 
  % so far, store it (and the matching points), else ignore
  if sum(curr) > best
	best = sum(curr);
	loc = curr;
  end
end
x1_curr = x1(loc,:);
x2_curr = x2(loc,:);

% Recompute the homography using ALL of the best set 
% of points found above
q = 1;
  for c = 1:size(x1_curr,1)
    x_h = x1_curr(c,1); 
    xprime_h = x2_curr(c,1); 
    y_h = x1_curr(c,2);
    yprime_h = x2_curr(c,2);
    
    A_curr(q,:) = [x_h y_h 1 0 0 0 -xprime_h*x_h -xprime_h*y_h -xprime_h];
    A_curr(q+1,:) = [0 0 0 x_h y_h 1 -yprime_h*x_h -yprime_h*y_h -yprime_h];

    q = q+2;
  end
  
  % must find A'A with smallest eigenvalue (need diagonalized matrix) to
  % get homography matrix, and reshape homography matrix to 3x3
  [~,~,H_curr] = svd(A_curr);
  H_curr = reshape(H_curr(:,end),3,3);
  %make homography projective
  %T = maketform('projective',H_curr);
  T = H_curr;