function [square_x, square_y] = GetCoordinates(image)
im = imread(image);
imshow(im);
[square_x, square_y] = ginput()