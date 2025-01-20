function PIVParams = readPIVParameters(ParametersFileName)

%
% readPIVParameters - read the specified parameters file and assign values
% to the structure PIVParams.
% (This function is called by the function PIV)
%

%% Defaults
PIVParams = PIV_setDefaults();

%% Status settings
% These are used by the PIV script to remember in which part of the process
% it is
% Are we currently in the grid refinement procedure:
PIVParams.InGridRefinement = false;


%% Read parameters file
fid = fopen(ParametersFileName);    %Open file
line = fgetl(fid);    %Read first line
while (isempty(line))    %If first line is empty, read next line, and continue until we find a non-empty line
    line = fgetl(fid);
end
% Evaluate first non-empty line:
if (~strncmp(line, '#PIV-Parameters', 11) )    %Check if it is the right file (first line contains '#Experiment')
    error('Parameters file does not start with #PIV-Parameters');
end
line = fgetl(fid);
while (ischar(line))
    if (~isempty(line))
        if (line(1)~='%')    %Neglect empty and commented lines
            [token, rem] = strtok(line, ' ');
            while (rem(1)==' ')
                rem = rem(2:end);       %remove leading space
                if numel(rem)==0    % prevent error in case the field is empty
                    break
                end
            end
            switch token
                
                
                
                %%%%%%%%%%%%%
                %%% Tasks %%%
                %%%%%%%%%%%%%
                %%% Tasks: what do you want the PIV-program to do?
                % Do PIV analysis (cross correlation and stuff) if you choose not to, the
                % program will load data from the position where it would otherwise write
                % it
                case '#DoPIVAnalysis:'
                    if rem(1)=='t'
                        PIVParams.DoPIVAnalysis = true;
                    else
                        PIVParams.DoPIVAnalysis = false;
                    end
                % Convert dx's and dy's to u's and v's in physical units
                case '#ConvertData:'
                    if rem(1)=='t'
                        PIVParams.ConvertData = true;
                    else
                        PIVParams.ConvertData = false;
                    end
                % Calculate scalar fields:
                case '#CalcScalars:'
                    if rem(1)=='t'
                        PIVParams.CalcScalars = true;
                    else
                        PIVParams.CalcScalars = false;
                    end
                % Use existing data files only for reading, otherwise they will be
                % overwritten
                case '#SaveData:'
                    PIVParams.SaveData = eval(rem);
                % Display the vectors of every frame that is processed or loaded
                % (last image is always displayed)
                case '#DisplayEveryFrame:'
                    PIVParams.DisplayEveryFrame = eval(rem);
                % Save images (only displayed images are saved)
                case '#SaveImages:'
                    PIVParams.SaveImages = eval(rem);
                
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% File locations and names %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % filename prefix, e.g. 'img0001.png' has ImageFileName 'img'
                case '#ImageFileName:'
                    PIVParams.ImageFileName = rem;
                % filename extension, e.g. 'img0001.png' has extension 'png'
                case '#FileExtension:'
                    PIVParams.FileExtension = rem;
                % are the images in a multipage tiff?
                case '#MultipageTiff:'
                    PIVParams.MultipageTiff = eval(rem);
                % bit depth of the images (default is 8)
                case '#ImageBitDepth:'
                    PIVParams.ImageBitDepth = str2double(rem);
                % The directory where the images are located:
                case '#Directory:'
                    if (rem(end) ~= '/')
                        rem = [rem '/']; %#ok<AGROW>
                    end
                    PIVParams.Directory = rem;
                % The directory where the PIV-output images will be placed:
                case '#ImageOutputDir:'
                    if (rem(end) ~= '/')
                        rem = [rem '/']; %#ok<AGROW>
                    end
                    PIVParams.ImageOutputDir = rem;
                % the file format for the output images
                case '#ImageOutputFormat:'
                    PIVParams.ImageOutputFormat = rem;
                % The directory where the PIV-output data will be saved or read:
                case '#DataDir:'
                    if (rem(end) ~= '/')
                        rem = [rem '/']; %#ok<AGROW>
                    end
                    PIVParams.DataDir = rem;
                % The first image that will be processed (a number).
                case '#FirstImage:'
                    PIVParams.FirstImage = str2double(rem);
                % The last image that will be processed
                case '#LastImage:'
                    PIVParams.LastImage = str2double(rem);
                % The number of images we want to skip between to PIV images (set to 0 (zero) if you want to process all images)
                case '#SkipImages:'
                    PIVParams.SkipImages = str2double(rem);
                    
                    
                    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% pre-processing settings %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                case '#UseMinMax:'
                    PIVParams.UseMinMax = eval(rem);
                case '#MinMaxKernel:'
                    PIVParams.MinMaxKernel = str2double(rem);
                case '#UseUnsharpMask:'
                    PIVParams.UseUnsharpMask = eval(rem);
                case '#UnsharpMaskAmount:'
                    PIVParams.UnsharpMaskAmount = str2double(rem);
                case '#UseGaussianBlur:'
                    PIVParams.UseGaussianBlur = eval(rem);
                case '#GaussianBlurRadius:'
                    PIVParams.GaussianBlurRadius = str2double(rem);
                    
                    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% ROI and grid settings %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % rotate the image counterclockwise by x degrees (default = 0)
                case '#RotationAngle:'
                    PIVParams.RotationAngle = str2double(rem);
                %%% Use ROI to crop image or only to limit the region for
                %%% cross-correlations:
                case '#CropImageToROI:'
                    PIVParams.CropImageToROI = eval(rem);
                %%% static ROI
                % Select a static region of interest (ROI)
                case '#UseROI:'
                    PIVParams.UseROI = eval(rem);
                case '#ROILeftEdge:'
                    PIVParams.ROILeftEdge = str2double(rem);
                case '#ROIRightEdge:'
                    PIVParams.ROIRightEdge = str2double(rem);
                case '#ROITopEdge:'
                    PIVParams.ROITopEdge = str2double(rem);
                case '#ROIBottomEdge:'
                    PIVParams.ROIBottomEdge = str2double(rem);
                %%% moving ROI
                % Select a moving ROI
                case '#UseMovingROI:'
                    PIVParams.UseMovingROI = eval(rem);
                % initial left and right edge
                case '#MovingROILeftStart:'
                    PIVParams.MovingROILeftStart = str2double(rem);
                case '#MovingROIRightStart:'
                    PIVParams.MovingROIRightStart = str2double(rem);
                % horizontal velocity in pixels/frame (doen't have to be integer)
                case '#MovingROIHorizontalVelocity:'
                    PIVParams.MovingROIHorizontalVelocity = str2double(rem);
                % initial top and bottom edge
                case '#MovingROITopStart:'
                    PIVParams.MovingROITopStart = str2double(rem);
                case '#MovingROIBottomStart:'
                    PIVParams.MovingROIBottomStart = str2double(rem);
                % vertical velocity in pixels/frame (doen't have to be integer)
                case '#MovingROIVerticalVelocity:'
                    PIVParams.MovingROIVerticalVelocity = str2double(rem);
                %%% Grid size and overlap
                % correlation grid size (size of correlation windows, in pixels)
                case '#GridSize:'
                    PIVParams.GridSize = str2double(rem);
                % overlap of the correlation windows (in pixels)
                case '#Overlap:'
                    PIVParams.Overlap = str2double(rem);
                
                     
                    
                    
                    
                    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Correlation settings %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % method of correlation: XCORR or NORMXCORR
                % 'XCORR' is faster, but the correlation coefficient does not have values
                % between -1 and 1, which makes it difficult to use as a quality indicator
                % 'NORMXCORR' is slower (the few tests I did were 2x slower), but having a
                % nice quality indicator may be well worth it. I recommend using this one,
                % and then you can for example say that you do not trust vectors with a 
                % correlation coefficient less than 0.5.
                case '#Method:'
                    PIVParams.Method = rem;
                % number of passes (not for gridrefinement)
                case '#NumPasses:'
                    PIVParams.NumPasses = str2double(rem);
                % use grid refinement scheme
                case '#UseGridRefinement:'
                    PIVParams.UseGridRefinement = eval(rem);
                % Initial grid size for grid refinement (should be power of 2)
                case '#InitialGridSize:'
                    PIVParams.InitialGridSize = str2double(rem);
                % Correlation step: correlate image i with i+n (n larger or
                % equal to 1)
                case '#CorrelationStep:'
                    PIVParams.CorrelationStep = str2double(rem);

                    
                    
                    
                    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Post-processing settings %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % maximum number of validation passes (only makes sense if you also set 
                % some outlier replacement scheme)
                case '#numValidationPasses:'
                    PIVParams.numValidationPasses = str2double(rem);
                % use maximum displacement vector, if vector is larger, the vector is removed
                case '#UseMaxDisplacement:'
                    PIVParams.UseMaxDisplacement = eval(rem);
                % set the maximum displacement vector (maximum length of a
                % vector in pixels)
                case '#MaxDisplacement:'
                    PIVParams.MaxDisplacement = str2double(rem);
                % use minimum correlation coefficient (NORMXCORR method recommende, see
                % method above)
                case '#UseMinCorrelationCoefficient:'
                    PIVParams.UseMinCorrelationCoefficient = eval(rem);
                % set the minimum correlation coefficient
                case '#MinCorrelationCoefficient:'
                    PIVParams.MinCorrelationCoefficient = str2double(rem);
                % use the normalized median test (should be a very good test)
                case '#UseNormalizedMedianTest:'
                    PIVParams.UseNormalizedMedianTest = eval(rem);
                % set the threshold for the normalized median test
                case '#MedianThreshold:'
                    PIVParams.MedianThreshold = str2double(rem);
                %%% Outlier replacement
                % try second peak (this makes the code significantly slower, so test to see
                % if you really need it)
                case '#UseSecondPeak:'
                    PIVParams.UseSecondPeak = eval(rem);
                % interpolate spurious data
                case '#UseInterpolation:'
                    PIVParams.UseInterpolation = eval(rem);
                %%% Data smoothing
                % smooth data
                case '#UseDataSmoothing:'
                    PIVParams.UseDataSmoothing = eval(rem);
                % set gaussian smoothing filter size
                case '#SmoothingFilterSize:'
                    PIVParams.SmoothingFilterSize = str2double(rem);
                case '#SmoothingFilterType:'
                    PIVParams.SmoothingFilterType = rem;
                    
                
                    
                    
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Conversion factors %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                % Framerate (frames per second)
                case '#FPS:'
                    PIVParams.FPS = str2double(rem);
                % Scaling (meters per pixel)
                case '#Scaling:'
                    PIVParams.Scaling = str2double(rem);
                % Units
                case '#TimeUnit:'
                    PIVParams.TimeUnit = rem;
                case '#LengthUnit:'
                    PIVParams.LengthUnit = rem;
                case '#VelocityUnit:'
                    PIVParams.VelocityUnit = rem;
                

                    
                    
                    
                %%%%%%%%%%%%%%%%%%%%%%%%
                %%% Display settings %%%
                %%%%%%%%%%%%%%%%%%%%%%%%
                % Show the vectors:
                case '#ShowVectors:'
                    PIVParams.ShowVectors = eval(rem);
                % Color code the vectors:
                case '#ColorCodeVectors:'
                    PIVParams.ColorCodeVectors = eval(rem);
                case '#VectorColorValidation:'
                    PIVParams.VectorColorValidation = eval(rem);
                case '#VectorColormap:'
                    PIVParams.VectorColormap = rem;
                % Color of the vectors (yellow, blue, red, green, white, black)
                % Or in rgb format: [r g b] with r,g,b between 0 and 1
                case '#VectorColor:'
                    if rem(1)=='['
                        eval(['PIVParams.VectorColor = ' rem])
                    else
                        PIVParams.VectorColor = rem;
                    end
                % The scale of the displayed vectors (pixels per m/s)
                case '#VectorScale:'
                    PIVParams.VectorScale = str2double(rem);
                % Show the image as a background for the vectors
                case '#ShowImage:'
                    PIVParams.ShowImage = eval(rem);
                case '#ShowGrid:'
                    % Show the grid of correlation windows
                    PIVParams.ShowGrid = eval(rem);
                case '#ShowContour:'
                    % Show a scalar field as contour
                    PIVParams.ShowContour = eval(rem);
                case '#ContourScalar:'
                    % Which scalar? (magnitude, signalToNoiseRatio, correlationCoefficient,
                    % vorticity, shearStrain, normalStrain, xComponent, yComponent, dudx, dudy,
                    % dvdx, dvdy)
                    PIVParams.ContourScalar = rem;
                case '#ContourColormap:'
                    PIVParams.ContourColormap = rem;
                case '#ContourSteps:'
                    PIVParams.ContourSteps = str2double(rem);
                case '#DisplayScalarInfo:'
                    % Display global information sbout the scale (min, max, mean value)
                    PIVParams.DisplayScalarInfo = eval(rem);
                % Set the range of values to be displayed of the scalar field (or use autoscale)
                case '#ScalarMinValue:'
                    PIVParams.ScalarMinValue = str2double(rem);
                case '#ScalarMaxValue:'
                    PIVParams.ScalarMaxValue = str2double(rem);
                case '#ScalarAutoScale:'
                    PIVParams.ScalarAutoScale = eval(rem);
                case '#ScalarPrefactor:'
                    PIVParams.ScalarPrefactor = str2double(rem);
                case '#OutputImageWidth:'
                    PIVParams.OutputImageWidth = str2double(rem);
                case '#OutputImageHeight:'
                    PIVParams.OutputImageHeight = str2double(rem);
                case '#FontName:'
                    PIVParams.FontName = rem;
                case '#FontSize:'
                    PIVParams.FontSize = str2double(rem);
                case '#ScalarLabel:'
                    if numel(rem)==0
                        PIVParams.ScalarLabel = '';
                    else
                        PIVParams.ScalarLabel = rem;
                    end

                    
                    
                    
                %%%%%%%%%%%%%%%%%%%%%%%%
                % if the property is not recognized:
                otherwise
                    display(['WARNING: Invalid property: ' line]);
            end
        end
    end
    line = fgetl(fid);    %Done with this one, read the next line...
end
fclose(fid);
