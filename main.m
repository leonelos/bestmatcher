%% Loading the data
% Available tests: 'barrsmith', 'booksh', 'Kyoto', and 'johnssonb'
test = 'Kyoto';
data = dlmread(strcat('data/', test, '_kps.txt'));
databy7 = dlmread(strcat('data/inliers_by7.txt'));
databy8 = dlmread(strcat('data/inliers_by8.txt'));
databyRansac = dlmread(strcat('data/inliers_byRansac.txt'));
databyLmeds = dlmread(strcat('data/inliers_byLmeds.txt'));

Fby7 = dlmread(strcat('data/fundamentalM_by7.txt'));
Fby8 = dlmread(strcat('data/fundamentalM_by8.txt'));
FbyRansac = dlmread(strcat('data/fundamentalM_byRansac.txt'));
FbyLmeds = dlmread(strcat('data/fundamentalM_byLmeds.txt'));

n = size(data, 1);
maxvalues = max(data(:,1:2));

%% RANSAC
max_iterations = inf;
iterations = 2000;
threshold = 3.0;
best_inliers = [];
confidence = 0.99;
bestF = []; %the best fundamental matrix



for i = 1 : iterations
    indices = randperm(n, 5);
    pts1 = data(indices, 1:2);
    pts2 = data(indices, 3:4);
    alphas = data(indices, 5);
    
    % Estimate a homography from the first three correspondences
    [normedPts1, normedPts2, T1, T2] = NormalizePoints(pts1, pts2);
    H = GetHomographyFromSIFT(normedPts1(1:3, :), normedPts2(1:3, :), alphas(1:3, :));
    H = inv(T2) * H * T1;    
    
    % Check if the sstwo addition points lie on the same plane
    proj1 = H * [pts1(4:5,:), repmat(1, 2, 1)]';
    is_too_close = 0;
    for j = 1 : 2
        dist = norm(proj1(1:2,j) / proj1(3,j) - pts2(j + 3,:)');
        if dist < threshold
            is_too_close = 1;
            break;
        end
    end
    
    if is_too_close
        continue;
    end
    
    % Estimate a fundamental matrix from the homography and point correspondences
    F = GetFundamentalMatrixFromHomographies(H, pts1, pts2, maxvalues);
       
    if F == 0
        continue;
    end
    
    if size(F, 1) ~= 3
        continue;
    end    
    
    % Get the inliers
    for fi = 1 : size(F, 3)
        inliers = [];
        Fi = F(:,:,fi);
        % Estimate the symmetric epipolar distance for each correspondences
        for j = 1 : n 
            l1 = [data(j, 3:4), 1] * Fi;
            l2 = Fi * [data(j, 1:2), 1]';

            l1 = l1 / sqrt(l1(1)^2 + l1(2)^2);
            l2 = l2 / sqrt(l2(1)^2 + l2(2)^2);

            dist = abs(l1 * [data(j, 1:2), 1]') + abs([data(j, 3:4), 1] * l2) * 0.5;

            if dist < threshold
                inliers = [inliers j];
            end
        end
        
        if length(best_inliers) < length(inliers)
            % Update inliers of the so-far-the-best model
            best_inliers = inliers;
            
            % Update max iteration number
            max_iterations = log(1 - confidence) / log(1 - (length(best_inliers) / n)^5);
        
            %Here I need to do the refinement 
            best_A = []
            r = length(best_inliers); 
            for ii = 1 : r 
              m = data(best_inliers(ii),:); %match
              x1 = m(1);
              y1 = m(2);
              x2 = m(3);
              y2 = m(4);
              row = [x1*x2 y1*x2 x2 x1*y2 y1*y2 y2 x1 y1 1];
              best_A = [best_A; row];
            end
            %yeah I got it, now let's do the SVD from A... so easy
            [U, S, V] = svd(best_A,0);
            Vt = V';
            %Vector<float> fbest = Vt_Abest.getRow(8);
            fbest = Vt(9,:);
            F_best = [fbest(1:3); fbest(4:6); fbest(7:9)];
            [UF_best, SF_best, VF_best] = svd(F_best,0);
            
            %I force s3 to be 0
            SF_best_mat = [SF_best(1,1) 0 0; 0 SF_best(2,2) 0; 0 0 0];
            
            %and now lets recompose the fundamental matrix matrix 
            bestF = UF_best * SF_best_mat * VF_best'  
        end
    end
    
    if i > max_iterations
        break;
    end    
end

%% Visualization
fprintf('Number of loaded points = %d\n', n);
fprintf('Number of found inliers = %d\n', length(best_inliers));
fprintf('Number of iterations = %d\n', i);

close all;
img1 = imread(strcat('data/', test, 'A.jpg'));
img2 = imread(strcat('data/', test, 'B.jpg'));

image = [img1 img2];

imshow(image)
hold on;

inliers_data = [];

for i = 1 : length(best_inliers)
    inliers_data = [inliers_data; data(best_inliers(i),:)];
    color = rand(1,3);
    
    plot([data(best_inliers(i), 1), size(img1, 2) + data(best_inliers(i), 3)], [data(best_inliers(i), 2) data(best_inliers(i), 4)], 'Color', color)
    
    scatter(data(best_inliers(i), 1), data(best_inliers(i), 2), 40, color, 'filled');
    scatter(size(img1, 2) + data(best_inliers(i), 3), data(best_inliers(i), 4), 40, color,'filled');
end
hold off;

%Here I will calculate the fundamental matrix from all inliers :)

in_points1 = inliers_data(:, 1:2);
in_points2 = inliers_data(:, 3:4);

%I homogenize the points
in_points1H = [in_points1 ones(size(in_points1,1), 1)];
in_points2H = [in_points2 ones(size(in_points2,1), 1)];

%and then normalize xD
[in_normedPts1, in_normedPts2, in_T1, in_T2] = NormalizePoints(in_points1H, in_points2H)

%and now lets build the constraint matrix ;)
n_inliers = size(inliers_data,1);
x1 = in_normedPts1';
x2 = in_normedPts2';
A = [x2(1,:)'.*x1(1,:)'   x2(1,:)'.*x1(2,:)'  x2(1,:)' ...
     x2(2,:)'.*x1(1,:)'   x2(2,:)'.*x1(2,:)'  x2(2,:)' ...
     x1(1,:)'             x1(2,:)'            ones(n_inliers,1) ]; 
%let's calculate      
[U,D,V] = svd(A);  

% Extract fundamental matrix from the column of V corresponding to smallest singular value.
F = reshape(V(:,9),3,3)'

% Enforce constraint that fundamental matrix has rank 2 by performing
% a svd and then reconstructing with the two largest singular values.
[U,D,V] = svd(F,0);
F = U*diag([D(1,1) D(2,2) 0])*V';

% Denormalise
F = T2'*F*T1;

%F = bestF; %looks like my bestF doesnt work 

triangPointsby7 = triangulation(databy7, Fby7);
triangPointsby8 = triangulation(databy8, Fby8);
triangPointsbyRansac = triangulation(databyRansac, FbyRansac);
triangPointsbyLmeds = triangulation(databyLmeds, FbyLmeds);
triangPoints = triangulation(inliers_data, F);

figure('name','3D by 7'); plot3(triangPointsby7(:,1),triangPointsby7(:,2),-triangPointsby7(:,3)*100,'.');
figure('name','3D by 8'); plot3(triangPointsby8(:,1),triangPointsby8(:,2),-triangPointsby8(:,3)*100,'.');
figure('name','3D by Ransac'); plot3(triangPointsbyRansac(:,1),triangPointsbyRansac(:,2),-triangPointsbyRansac(:,3)*100,'.');
figure('name','3D by Lmeds'); plot3(triangPointsbyLmeds(:,1),triangPointsbyLmeds(:,2),-triangPointsbyLmeds(:,3)*100,'.');
figure('name','3D by 5 points'); plot3(triangPoints(:,1),triangPoints(:,2),-triangPoints(:,3)*100,'.');

%triangPoints
