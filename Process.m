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
    
    %Run Harris corner detector on two images and then Ransac on the two
    %corner matrices (default method is Harris_
    c1 = corner(I1);
    c2 = corner(I2);

    %Run RANSAC on corner matrices
    [corner_inliers1, corner_inliers2] = Ransac(c1, c2);