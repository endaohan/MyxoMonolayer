% *********************************************************************** %
% Read and write parameters for TFM PIV
% Written by Endao Han, version 5, modified 2020/5/25. 
% *********************************************************************** %

function FUN_PIV_readWriteParameters_TFM(route, Parameter, iRef, flagShow )

route_in = [pwd '\TFM_PIV']; 
if nargin == 3
    flagShow = 1; 
end


EX = exist(route,'file'); 
if EX == 7
    %% Get information of the images
    % Find the names of all images
    fileList = FUN_FileName([route '\Images'],'.tif','Img_','No'); 
    N_img = length(fileList); 
    % Find image size
    I = imread([route '\Images\' fileList{1,1} '.tif']); 
    dim = size(I); 
    
    %% Change the parameters
    % Read in the template
    PIVParams = PIV_readPIVParameters([route_in '\PIV_Parameters.txt']);
    
    PIVParams.DisplayEveryFrame = flagShow; 
    
    % Change the input and output ranges
    PIVParams.Directory = [route '\Images']; 
    PIVParams.ImageOutputDir = [route '\TFM_Images']; 
    PIVParams.DataDir = [route '\TFM_Data\TFM_Data_Ref' num2str(iRef)];
    
    % Number of frames
    if isfield(Parameter,'tRange')
        ind = find( Parameter.tRange(:,1) == iRef ); 
        PIVParams.FirstImage = Parameter.tRange(ind,2); 
        PIVParams.LastImage = Parameter.tRange(ind,3); 
    else
        PIVParams.FirstImage = 1; 
        PIVParams.LastImage = N_img; 
    end
    PIVParams.SkipImages = Parameter.PIVStep; 
    
    % Window
    PIVParams.ROILeftEdge = 21; 
    PIVParams.ROIRightEdge = dim(2)-20; 
    PIVParams.ROITopEdge = 21; 
    PIVParams.ROIBottomEdge = dim(1)-50; 
    % PIVParams.CropImageToROI = 'false'; 
    
    % PIVParams.CropImageToROI = 'true'; 
    
    % Grid size
    % PIVParams.GridSize = 40; 
    % PIVParams.Overlap = 20; 
    % PIVParams.GridSize = 30; 
    % PIVParams.Overlap = 20; 
    PIVParams.GridSize = 20; 
    PIVParams.Overlap = 10; 
    
    % Correlation settings
    PIVParams.NumPasses = 1; 
    
    
    % Conversion
    PIVParams.FPS = 1; 
    PIVParams.Scaling = Parameter.lscl(1)/1E6; 
    
    
    % Visualization
    PIVParams.VectorScale = 10/PIVParams.FPS/PIVParams.Scaling; 
    
    
    % Write out the parameters
    PIV_writePIVParameters(PIVParams, ... 
        [route '\parameters_TFM_Ref' num2str(iRef) '.txt']); 
    
    
end




