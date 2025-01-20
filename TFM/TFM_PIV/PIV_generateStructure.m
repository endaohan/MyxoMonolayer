function Vectors = PIV_generateStructure(PIVParams, Im1)

%
% PIV_generateStructure - create the fields of the matrices (filled with
% zeros) that will be filled with data later
% (This function is called by the function PIV)
%


%% First determine the number of vectors N that we will generate


step = PIVParams.GridSize - PIVParams.Overlap;
grid = PIVParams.GridSize; %This is a bit shorter than everytime PIVParams.GridSize :)

% get the size of the images
[height, width] = size(Im1);

if PIVParams.UseROI
    % we have a helper function which creates the correct coordinates for us:
    [x_points, y_points, ~] = PIV_ROI(PIVParams);
%     % set points only inside the Region Of Interest
%     y_points = PIVParams.ROITopEdge:step:PIVParams.ROIBottomEdge-grid;
%     x_points = PIVParams.ROILeftEdge:step:PIVParams.ROIRightEdge-grid;
else
    % all the corners of our correlation windows
    y_points = 1:step:height-grid+1;
    x_points = 1:step:width-grid+1;
end

% calculate the total number of coorelation windows
N = length(y_points) * length(x_points);

%% pre-allocate matrices
% This is the main purpose of this function
Vectors.x = zeros(1,N);
Vectors.y = zeros(1,N);
Vectors.dx = zeros(1,N);
Vectors.dy = zeros(1,N);
Vectors.dx2 = zeros(1,N);
Vectors.dy2 = zeros(1,N);
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
