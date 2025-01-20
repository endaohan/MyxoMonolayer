% *********************************************************************** %
% Read and write parameters for TFM PIV, for the slow drift. 
% Written by Endao Han, Version 4, modified 2020/5/25. 
% *********************************************************************** %

function FUN_PIV_readWriteParameters_LowRes(route, Parameter)

% Laptop (0) or lab desktop (1)
flag_computer = 1; 
if flag_computer == 0
    route_in = 'E:\Projects\2018_FruitingBody\Code'; 
elseif flag_computer == 1
    route_in = 'E:\OneDrive - Princeton University\Projects\2018_FruitingBody\Code'; 
end


EX = exist(route,'file'); 
if EX == 7
    %% Get information of the images
    % Find the names of all images
    fileList = FUN_FileName([route '\Images_LR'],'.tif','Img_','No'); 
    N_img = length(fileList); 
    % Find image size
    % I = imread([route '\Images\' fileList{1,1} '.tif']); 
    % dim = size(I); 
    
    %% Change the parameters
    % Read in the template
    PIVParams = PIV_readPIVParameters([route_in '\PIV_Parameters.txt']);
    
    % Change the input and output ranges
    PIVParams.Directory = [route '\Images_LR'];
    PIVParams.ImageOutputDir = [route '\LowRes_Images'];
    PIVParams.DataDir = [route '\LowRes_Data'];
    
    % Number of frames    
    if isfield(Parameter,'tRange')
        PIVParams.FirstImage = Parameter.tRange(1,2); 
        PIVParams.LastImage = Parameter.tRange(end,3); 
    else
        PIVParams.FirstImage = 1; 
        PIVParams.LastImage = N_img; 
    end
    PIVParams.SkipImages = Parameter.PIVStep; 
    
    % Window
    PIVParams.ROILeftEdge = Parameter.bd_LM(1,1); 
    PIVParams.ROIRightEdge = Parameter.bd_LM(2,1); 
    PIVParams.ROITopEdge = Parameter.bd_LM(1,2); 
    PIVParams.ROIBottomEdge = Parameter.bd_LM(2,2); 
    % PIVParams.CropImageToROI = 'false'; 
    
    % PIVParams.CropImageToROI = 'true'; 
    
    % Grid size
    PIVParams.GridSize = 400; 
    PIVParams.Overlap = 200; 
    
    % Conversion
    PIVParams.FPS = 1; 
    % PIVParams.Scaling = Parameter.lscl(1)/1E6; 
    PIVParams.Scaling = 1; 
    
    
    % Visualization
    PIVParams.VectorScale = 5/PIVParams.FPS/PIVParams.Scaling; 
    
    
    % Write out the parameters
    PIV_writePIVParameters(PIVParams, [route '\parameters_LowRes.txt']); 
    
    
end




