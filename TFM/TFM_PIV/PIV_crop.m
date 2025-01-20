function [PIVParams, Im] = PIV_crop(PIVParams, Im)

%
% PIV_crop - cropping function for PIV
% (This function is called by the function PIV and PIV_preprocess)
%


%% Crop image
if PIVParams.CropImageToROI
    % in case we're using a moving ROI, we have to update PIVParams for the
    % current ROI position:
    [~, ~, PIVParams] = PIV_ROI(PIVParams);
   
    % crop the image
    if PIVParams.UseROI
        Top = PIVParams.ROITopEdge;
        Bottom = PIVParams.ROIBottomEdge;
        Left = PIVParams.ROILeftEdge;
        Right = PIVParams.ROIRightEdge;
        [h, w] = size(Im);
        % Check if we can actually crop the image:
        if Bottom>h
            PIVParams.Abort = true;
            PIVParams.AbortMessage = 'PIV_proprocess.m Code E1 - ROI outside image';
            PIVParams.ROIBottomEdge = h;
            Bottom = h;
        end
        if Right>w
            PIVParams.Abort = true;
            PIVParams.AbortMessage = 'PIV_proprocess.m Code E2 - ROI outside image';
            PIVParams.ROIRightEdge = w;
            Right = w;
        end
        % Do the crop:
        Im = Im(Top:Bottom, Left:Right);
    end
end