function [x_points, y_points, PIVParams] = PIV_ROI(PIVParams)

% PIV_ROI
% [x_points, y_points] = PIV_ROI(PIVParams)
%
% makes an array of x_points and y_points according to the ROI (Region Of
% Interest)
% This can be a static or a moving (position of the ROI depends on the
% of the image) ROI.
%

%% adjust current ROI if the ROI is moving:
if PIVParams.UseMovingROI
    % get time (= image number)
    t = PIVParams.imageNumber;
    
    % initial horizontal position:
    leftEdge = PIVParams.MovingROILeftStart;
    rightEdge = PIVParams.MovingROIRightStart;
    % horizontal displacement:
    displacement = round(t * PIVParams.MovingROIHorizontalVelocity);
    % current position
    PIVParams.ROILeftEdge = leftEdge + displacement;
    PIVParams.ROIRightEdge = rightEdge + displacement;
    
    % initial vertical position:
    topEdge = PIVParams.MovingROITopStart;
    bottomEdge = PIVParams.MovingROIBottomStart;
    % vertical displacement:
    displacement = round(t * PIVParams.MovingROIVerticalVelocity);
    % current position
    PIVParams.ROITopEdge = topEdge + displacement;
    PIVParams.ROIBottomEdge = bottomEdge + displacement;
end

%% set x_points and y_points
% grid size:
grid = PIVParams.GridSize;
% step size:
step = PIVParams.GridSize - PIVParams.Overlap;



if PIVParams.CropImageToROI
    width = 1 + PIVParams.ROIRightEdge - PIVParams.ROILeftEdge;
    height = 1 + PIVParams.ROIBottomEdge - PIVParams.ROITopEdge;
    x_points = 1:step:width-grid+1;
    y_points = 1:step:height-grid+1;
else
    % set points inside the ROI (Region Of Interest)
    x_points = PIVParams.ROILeftEdge:(PIVParams.GridSize - PIVParams.Overlap):PIVParams.ROIRightEdge-PIVParams.GridSize;
    y_points = PIVParams.ROITopEdge:(PIVParams.GridSize - PIVParams.Overlap):PIVParams.ROIBottomEdge-PIVParams.GridSize;
end

