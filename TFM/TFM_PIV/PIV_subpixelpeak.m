function [deltaY, deltaX] = PIV_subpixelpeak(I)

%
% PIV_subpixelpeak - Find the subpixel position of a correlation peak
% (This function is called by the function PIV)
%

% make sure we have no negative values
I = I - min(I(:));

% also, we don't want to be taking the log of exactly zero. Simple fix:
I = I + 1;

% calculate the subpixel peak position with a 3-point gaussian fit:
I = log(I);
deltaY = (I(1,2)-I(3,2))/(2*I(1,2)-4*I(2,2)+2*I(3,2));
deltaX = (I(2,1)-I(2,3))/(2*I(2,1)-4*I(2,2)+2*I(2,3));