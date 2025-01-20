function [PIVParams, Vectors] = PIV_refineGrid(PIVParams, Vectors, Im1, NewGridSize)

%
% PIV_refineGrid - double the density of correlation windows. The already
% existing vectors are interpolated to the new grid.
% (This function is called by the function PIV)
%


% read and copy vector data
x = Vectors.x;
y = Vectors.y;
dx = Vectors.dx;
dy = Vectors.dy;

% get dimensions of vector field
i = find(y>y(1),1,'first');
width = i-1;
height = length(x) / width;

% reshape vectors into matrices:
x = reshape(x,width,height)';
y = reshape(y,width,height)';
dx = reshape(dx,width,height)';
dy = reshape(dy,width,height)';

%% generate refined grid:
PIVParams.GridSize = NewGridSize;

step = PIVParams.GridSize - PIVParams.Overlap;
grid = PIVParams.GridSize; %This is a bit shorter than everytime PIVParams.GridSize :)

% get the size of the images
[height, width] = size(Im1);

% set the positions of the corners of our correlation windows:
if PIVParams.UseROI
    % we have a helper function which creates the correct coordinates for us:
    [x_points, y_points, ~] = PIV_ROI(PIVParams);
else
    % use the complete image
    y_points = 1:step:height-grid+1;
    x_points = 1:step:width-grid+1;
end

%% The loop going through all the correlation windows
index = 1;
for yy = y_points
    for xx = x_points
        % register current position:
        Vectors.y(index) = yy + round(grid/2);
        Vectors.x(index) = xx + round(grid/2);
        % increase index
        index = index + 1;
    end
end
%% Read in our refined grid

% read and copy vector data
xx = Vectors.x;
yy = Vectors.y;

% get dimensions of vector field
i = find(yy>yy(1),1,'first');
width = i-1;
height = length(xx) / width;

% reshape vectors into matrices:
xx = reshape(xx,width,height)';
yy = reshape(yy,width,height)';


%% interpolate existing data
dx = interp2(x,y,dx,xx,yy,'spline');
dy = interp2(x,y,dy,xx,yy,'spline');


% reshape matrices to vectors to be compatible to the other data
% number of vectors N:
N = width * height;
Vectors.x = reshape(xx',1,N);
Vectors.y = reshape(yy',1,N);
Vectors.dx = reshape(dx',1,N);
Vectors.dy = reshape(dy',1,N);
% the other fields are not important, but must have the same size, so we
% simply fill them with zero's. They will get their data later
Vectors.u = zeros(1,N);
Vectors.v = zeros(1,N);
Vectors.magnitude = zeros(1,N);
Vectors.signalToNoiseRatio = zeros(1,N);
Vectors.correlationCoefficient = zeros(1,N);
Vectors.vorticity = zeros(1,N);
Vectors.shearStrain = zeros(1,N);
Vectors.normalStrain = zeros(1,N);
Vectors.dudx = zeros(1,N);
Vectors.dudy = zeros(1,N);
Vectors.dvdx = zeros(1,N);
Vectors.dvdy = zeros(1,N);
Vectors.validationFlag = zeros(1,N);
