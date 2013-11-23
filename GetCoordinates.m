function [square_x, square_y] = GetCoordinates(image)
im = imread(image);
imshow(im);
format long g;
[square_x, square_y] = ginput();