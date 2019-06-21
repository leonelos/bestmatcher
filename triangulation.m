function [triangPoints] = triangulation(inliers_data, F)
      
  [U,D,V] = svd(F,0);
  e1 = hnormalise(V(:,3));
  e2 = hnormalise(U(:,3));

  %Now I calculate the skew symetrix matrix from e2
  v = e2;
  e_x = [ 0   -v(3)  v(2)
      v(3)   0   -v(1)
     -v(2)  v(1)   0 ];
     
  %Now! Let's calculate the projection matrix and then K and rotations from projection :)
  [U,S,V] = svd(F); %svd(F_best);
  e = U(:,3);
  P = [-vgg_contreps(e)*F e] %[-vgg_contreps(e)*F_best e]
  [K, Rc_w, Pc, pp, pv] = decomposecamera(P);
  %P = K*[Rc_w -Rc_w*Pc] %yeah it works

  %Let's do the triangulation
  x1 = cart_2_homo(inliers_data(:,1:2));
  x2 = cart_2_homo(inliers_data(:,3:4));
  numMatches = size(x1,1);
  triangPoints = zeros(numMatches, 3);
  projPointsImg1 = zeros(numMatches, 2);
  projPointsImg2 = zeros(numMatches, 2);

  %calcualte the triangulated points, + their projections onto each img plane
  camMatrix1 = eye(3,4)
  camMatrix2 = P
  camCenter1 = get_cam_center(camMatrix1);
  camCenter2 = get_cam_center(camMatrix2);
  for i = 1:numMatches
    pt1 = x1(i,:);
    pt2 = x2(i,:);
    crossProductMat1 = [  0   -pt1(3)  pt1(2); pt1(3)   0   -pt1(1); -pt1(2)  pt1(1)   0  ];
    crossProductMat2 = [  0   -pt2(3)  pt2(2); pt2(3)   0   -pt2(1); -pt2(2)  pt2(1)   0  ];    
    Eqns = [ crossProductMat1*camMatrix1; crossProductMat2*camMatrix2 ];
    
    [~,~,V] = svd(Eqns);
    triangPointHomo = V(:,end)'; %4 dim (3 dimensions + homo coord)
    %save the triangulated 3d point
    triangPoints(i,:) = homo_2_cart(triangPointHomo);
    
    %project the triangulated point using both camera matrices for later residual calculations
    projPointsImg1(i,:) = homo_2_cart((camMatrix1 * triangPointHomo')');
    projPointsImg2(i,:) = homo_2_cart((camMatrix2 * triangPointHomo')');
  end

  % plot the triangulated points and the camera centers
  %%plot_triangulation(triangPoints, camCenter1, camCenter2);

  %%matches = inliers_data;

  %calculate the error distance between the triangulated point projected onto
  %the image plane and the actual location of the point on the image plane
  %%distances1 = diag(dist2(matches(:,1:2), projPointsImg1));
  %%distances2 = diag(dist2(matches(:,3:4), projPointsImg2));
  %%display(['Mean Residual 1: ', num2str(mean(distances1))]);
  %%display(['Mean Residual 2: ', num2str(mean(distances2))]);
  
end
