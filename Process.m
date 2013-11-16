function [x1, x2] = Process(im1, im2)
% Takes two images as input, runs sift and then RANSAC
    % Compute and match SIFT descriptors (modified from sift_demo.m)
    %im1 = double(im1); im2 = double(im1);
    im1 = imread(im1); im2 = imread(im2);
    im1 = rgb2gray(im1); im2 = rgb2gray(im2);
    im1 = double(im1)/246; im2 = double(im2)/246;
    I1 = im1; I2 = im2;
    [frames1,descr1,gss1,dogss1] = sift( I1, 'Verbosity', 1 ) ;
    [frames2,descr2,gss2,dogss2] = sift( I2, 'Verbosity', 1 ) ;

    descr1=uint8(512*descr1) ;
    descr2=uint8(512*descr2) ;
    matches=siftmatch( descr1, descr2 ) ;

    figure(3) ; clf ;
    plotmatches(I1,I2,frames1(1:2,:),frames2(1:2,:),matches) ;
    drawnow ;

    x1 = frames1(1:2,matches(1,:));
    x2 = frames2(1:2,matches(2,:));
    
    x1 = [x1(1,:)', x1(2,:)'];
    x2 = [x2(1,:)', x2(2,:)'];
    
    %Run Harris corner detector on two images and then Ransac on the two
    %corner matrices (default method is Harris_
    c1 = corner(I1);
    c2 = corner(I2);

    %Run RANSAC on corner matrices
    [corner_inliers1, corner_inliers2] = Ransac(c1, c2);
    
    x1 = [x1; corner_inliers1];
    x2 = [x2; corner_inliers2];
    
    % Run RANSAC
    %[sift_r1, sift_r2, H] = Ransac(x1, x2);
    
    % Estimate fundamental matrix
    %F = estimateFundamentalMatrix(sift_r1,sift_r2);
    F = estimateFundamentalMatrix(x1,x2,'Method', 'RANSAC', 'NumTrials', 200, 'DistanceThreshold', 10); 
    sift_r1 = x1;
    sift_r2 = x2;
    %Set intrinsic camera matrix
    %K = [1138.81, 0, 535.107; 0, 1159.81, 298.384; 0, 0, 1];
    K = [832.85, 0.1401, 304.18; 0, 832.90, 206.76; 0, 0, 1];
    
    
    %Essential Matrix
    E = K'*F*K;
    [U,S,V] = svd(E);
    S_prime = zeros(3,3);
    
    
    diagonals = [S(1,1), S(2,2), S(3,3)];
    min_index = -1;
    for i=1:3
        if min_index == -1
            min_index = i;
        else
            if diagonals(i) < diagonals(min_index)
                min_index = i;
            end
        end
    end
    if min_index == 1
        avg = (diagonals(2) + diagonals(3))/2;
        S_prime(2, 2) = avg;
        S_prime(3, 3) = avg;
    elseif min_index == 2
        avg = (diagonals(1) + diagonals(3))/2;
        S_prime(1, 1) = avg;
        S_prime(3, 3) = avg;
    else
        avg = (diagonals(1) + diagonals(2))/2;
        S_prime(1, 1) = avg;
        S_prime(2, 2) = avg;
    end
    
    
    E_prime = U*S_prime*V';
    
    [U, S, V] = svd(E_prime);

    W = [0, -1, 0; 1, 0, 0; 0, 0, 1];
    
    %Possible R/t for second image project matrix (P2)
    u3 = U(:, 3);
    first = K * [U*W*V', u3];
    second = K * [U*W*V', -u3];
    third = K * [U*W'*V', u3];
    fourth = K * [U*W'*V', -u3];
    
    P1 = [K, zeros(3,1)];    

    scene_points_first = [];
    scene_points_second = [];
    scene_points_third = [];
    scene_points_fourth = [];
    for j = 1:size(sift_r1,1)
        A_first = [P1(3,:) * sift_r1(j,1) - P1(1,:); P1(3,:) * sift_r1(j,2) - P1(2,:); first(3,:) * sift_r2(j,1) - first(1,:); first(3,:) * sift_r2(j,2) - first(2,:)];
        A_second = [P1(3,:) * sift_r1(j,1) - P1(1,:); P1(3,:) * sift_r1(j,2) - P1(2,:); second(3,:) * sift_r2(j,1) - second(1,:); second(3,:) * sift_r2(j,2) - second(2,:)];
        A_third = [P1(3,:) * sift_r1(j,1) - P1(1,:); P1(3,:) * sift_r1(j,2) - P1(2,:); third(3,:) * sift_r2(j,1) - third(1,:); third(3,:) * sift_r2(j,2) - third(2,:)];
        A_fourth = [P1(3,:) * sift_r1(j,1) - P1(1,:); P1(3,:) * sift_r1(j,2) - P1(2,:); fourth(3,:) * sift_r2(j,1) - fourth(1,:); fourth(3,:) * sift_r2(j,2) - fourth(2,:)];
        
        [~,~,V_first] = svd(A_first);
        [~,~,V_second] = svd(A_second);
        [~,~,V_third] = svd(A_third);
        [~,~,V_fourth] = svd(A_fourth);
    
        % get last column of Vs and normalize
        last_col_first = V_first(:,end)/V_first(end,end);
        last_col_second = V_second(:,end)/V_second(end,end);
        last_col_third = V_third(:,end)/V_third(end,end);
        last_col_fourth = V_fourth(:,end)/V_fourth(end,end);
        
        scene_points_first = [scene_points_first, last_col_first];
        scene_points_second = [scene_points_second, last_col_second];
        scene_points_third = [scene_points_third, last_col_third];
        scene_points_fourth = [scene_points_fourth, last_col_fourth];
    end
    
    scene_points_first = [scene_points_first(1,:)',scene_points_first(2,:)',scene_points_first(3,:)',scene_points_first(4,:)'];
    scene_points_second = [scene_points_second(1,:)',scene_points_second(2,:)',scene_points_second(3,:)',scene_points_second(4,:)'];
    scene_points_third = [scene_points_third(1,:)',scene_points_third(2,:)',scene_points_third(3,:)',scene_points_third(4,:)'];
    scene_points_fourth = [scene_points_fourth(1,:)',scene_points_fourth(2,:)',scene_points_fourth(3,:)',scene_points_fourth(4,:)'];
    
    [x_grid, y_grid] = meshgrid(1:size(im1,2), 1:size(im1,1));
    
    scene_points_first = scene_points_first * 246;
    scene_points_second = scene_points_second * 246;
    scene_points_third = scene_points_third * 246;
    scene_points_fourth = scene_points_fourth * 246;
    
    num_neg_first = sum(scene_points_first(:,3) < 0);
    num_neg_second = sum(scene_points_second(:,3) < 0);
    num_neg_third = sum(scene_points_third(:,3) < 0);
    num_neg_fourth = sum(scene_points_fourth(:,3) < 0);
    
    if ( num_neg_first <= num_neg_second )
        if ( num_neg_third <= num_neg_fourth )
            scene_points_1 = scene_points_first;
            scene_points_2 = scene_points_third;
        else
            scene_points_1 = scene_points_first;
            scene_points_2 = scene_points_fourth;
        end
    elseif ( num_neg_second <= num_neg_first )
        if ( num_neg_third <= num_neg_fourth )
            scene_points_1 = scene_points_second;
            scene_points_2 = scene_points_third;
        else
            scene_points_1 = scene_points_second;
            scene_points_2 = scene_points_fourth;
        end
    end
    
    range_first = sqrt((max(scene_points_1(:,1)) - min(scene_points_1(:,1))^2) + (max(scene_points_1(:,2)) - min(scene_points_1(:,2))^2));
    range_second = sqrt((max(scene_points_2(:,1)) - min(scene_points_2(:,1))^2) + (max(scene_points_2(:,2)) - min(scene_points_2(:,2))^2));
    
    if ( range_first > range_second )
        scene_points = scene_points_1;
    else
        scene_points = scene_points_2;
    end
    
    z_i = griddata(scene_points(:,1), scene_points(:,2), scene_points(:,3), x_grid, y_grid, 'cubic');
    figure;
    warp(x_grid, y_grid, z_i, im1(:,:,1));
    view([0,90]);