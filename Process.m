function [x1, x2] = Process(im1, im2)
% Takes two images as input, runs sift and then RANSAC
    % Compute and match SIFT descriptors (modified from sift_demo.m)
    %im1 = double(im1); im2 = double(im1);
    im1 = imread(im1); im2 = imread(im2);
    I1 = rgb2gray(im1); I2 = rgb2gray(im2);
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
    
    % Run RANSAC
    [sift_r1, sift_r2, H] = Ransac(x1, x2);
    
    % Estimate fundamental matrix
    F = estimateFundamentalMatrix(sift_r1,sift_r2);
    Ftest = estimateFundamentalMatrix(x1,x2,'Method', 'RANSAC', 'NumTrials', 2000, 'DistanceThreshold', 0.49) 
    
    %Set intrinsic camera matrix
    K = [1138.81, 0, 535.107; 0, 1159.81, 298.384; 0, 0, 1];
    
    
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
    first = K * horzcat(U*W*V', u3);
    second = K * horzcat(U*W*V', -u3);
    third = K * horzcat(U*W'*V', u3);
    fourth = K * horzcat(U*W'*V', -u3);
    
    P1 = horzcat(K, zeros(3,1));
    
    %Run Harris corner detector on two images and then Ransac on the two
    %corner matrices (default method is Harris_
    c1 = corner(I1);
    c2 = corner(I2);

    %Run RANSAC on corner matrices
    [corner_inliers1, corner_inliers2] = Ransac(c1, c2);