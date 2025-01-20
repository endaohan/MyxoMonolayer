function TFM_PIV_ref(ParametersFileName,indRef)
%
% PIV - Perform PIV analysis based on a parameters file
% This is the only m-file that a user has to run
%
% usage:
% vectors = PIV(ParametersFileName)
%
% example:
% vectors = PIV('parameters.txt')
%
% PIV-software written by Ivo Peters, 2010-2014
% Modified by Endao Han for TFM data analysis on 2020/5/25
% Add one input parameter: indRef, which indicates the reference frame in
% the folder with all the images. 
%
% For debugging
% ParametersFileName = 'F:\DATA_Confocal\20200418\ExampleMove\parameters_TFM_Ref21.txt'; 

%% Set parameter
if (nargin == 1) || isempty(indRef)
    indRef = 1; 
end

% Display in command window? 1 - Yes, 0 - No. 
flag_disp = 0; 


%% Stuff we always do, to know what we're dealing with
disp('running PIV ...')

% make PIV parameters
PIVParams = readPIVParameters(ParametersFileName);

% make a parameter that we can use to abort PIV in case of a situation that
% would otherwise create an error. Anywhere in the program this parameter
% (PIVParams.Abort) can be set to true including a message with the reason 
% to abort. The message should at least have the name of the m-file where 
% the error occured and an identifier which makes it able for the user to 
% figure out exactly where the message was from. The message will be saved 
% to a file in the data output folder. 
% A typical reason to abort would be for example that the ROI is outside of
% the image.
PIVParams.Abort = false;
PIVParams.AbortMessage = '';

% check if the output directories exists, if not make them
if ~exist(PIVParams.DataDir,'dir') && PIVParams.SaveData
    mkdir(PIVParams.DataDir)
end
if ~exist(PIVParams.ImageOutputDir,'dir') && PIVParams.SaveImages
    mkdir(PIVParams.ImageOutputDir)
end

% make a list of the image files
PIVParams.files = dir([PIVParams.Directory PIVParams.ImageFileName '*.' PIVParams.FileExtension]);

% make a list of the data files (if they exist)
PIVParams.InputFiles = dir([PIVParams.DataDir 'PIV_VectorData' '*.mat']);

% make a figure handle so we know if a figure has been made and we can use
% it to refer to the figure
PIVParams.h_fig = 0;

tic
%% The big loop, going through all the images
for i = PIVParams.FirstImage:PIVParams.SkipImages+1:(PIVParams.LastImage-PIVParams.CorrelationStep)
    
    
    %% General stuff before every frame
    % Give the user some feedback:
    if flag_disp == 1
        disp('******************************')
        display(['Current frame: ' num2str(i)])
    end
    
    % save image number:
    PIVParams.imageNumber = i;
    
    % Data filename for current frame:
    DataFileName = [PIVParams.DataDir 'PIV_VectorData' num2str(i,'%06u') '.mat'];
        
    
    if PIVParams.MultipageTiff
        % read two images from multipage tiff file
        Im1 = imread([PIVParams.Directory PIVParams.files(1).name], i);
        Im2 = imread([PIVParams.Directory PIVParams.files(1).name], i+PIVParams.CorrelationStep);
    else
        % read two images from two subsequent files
        % For TFM, compare every frame with the first frame
        Im1 = imread([PIVParams.Directory PIVParams.files(indRef).name]);
        Im2 = imread([PIVParams.Directory PIVParams.files(i-1+PIVParams.CorrelationStep).name]);
    end
    
    % custom processing
    if exist('customProcess.m','file')
        Im1 = customProcess(Im1);
        Im2 = customProcess(Im2);
    end
    % rotate image
    Im1 = imrotate(Im1,PIVParams.RotationAngle);
    Im2 = imrotate(Im2,PIVParams.RotationAngle);
    
    % save original image for displaying later
    [PIVParams, ImOriginal] = PIV_crop(PIVParams, Im2);
    ImOriginal = ImOriginal/(2^(PIVParams.ImageBitDepth-8));
    
    % preprocess images
    [PIVParams, Im1, Im2] = PIV_prepocess_TFM(PIVParams, Im1, Im2);
    if PIVParams.Abort, break, end
    
    % pre-allocate matrices
    Vectors = PIV_generateStructure(PIVParams, Im1);
    
    %%
    if PIVParams.DoPIVAnalysis
        %% PIV Analysis
        %%% Grid refinement
        if PIVParams.UseGridRefinement
            [PIVParams, Vectors] = PIV_gridRefinement(PIVParams, Im1, Im2);
        end
                
        %%% Do the multipass correlation step
        %% override some user settings that we never/always use in multipass correlation
        % save user settings to recover after multipass step:
        %Overlap = PIVParams.Overlap;
        UseInterpolation = PIVParams.UseInterpolation;
        UseDataSmoothing = PIVParams.UseDataSmoothing;
        SmoothingFilterType = PIVParams.SmoothingFilterType;
        % and set the values we use for grid refinement:
        %PIVParams.Overlap = 0;
        PIVParams.UseDataSmoothing = true;
        PIVParams.SmoothingFilterType = 'average';
        PIVParams.UseInterpolation = true;
        %%% do the correlation and validation passes on target grid size
        if flag_disp == 1
            display(['Grid: ' num2str(PIVParams.GridSize)])
        end
        passNumber = 1;
        while passNumber<PIVParams.NumPasses
            dxOld = Vectors.dx;
            dyOld = Vectors.dy;
            display(['Pass: ' num2str(passNumber)])
            % correlate images
            Vectors = PIV_correlate(PIVParams, Vectors, Im1, Im2);
            % postprocess Vectors
            Vectors = PIV_postprocess(PIVParams, Vectors);
            % calculate change in dx's and dy's
            dxChange = max(abs(dxOld - Vectors.dx));
            dyChange = max(abs(dyOld - Vectors.dy));
            maxChange = max(dxChange, dyChange); 
            if flag_disp == 1
                display(['maxChange: ' num2str(maxChange)])
            end
            % if change is less than a pixel, we assume to have converged and
            % we can stop doing more passes
            if maxChange<1
                passNumber = PIVParams.NumPasses;
            end
            passNumber = passNumber + 1;
        end
        
        %%% last pass:
        %%% set back the replacement and filter settings:
        %PIVParams.Overlap = Overlap;
        PIVParams.UseInterpolation = UseInterpolation;
        PIVParams.UseDataSmoothing = UseDataSmoothing;
        PIVParams.SmoothingFilterType = SmoothingFilterType;
        
        dxOld = Vectors.dx;
        dyOld = Vectors.dy;
        if flag_disp == 1
            disp('Final pass')
        end
        % reset validation flag
        Vectors.validationFlag = 0*Vectors.validationFlag;
        % correlate images
        Vectors = PIV_correlate(PIVParams, Vectors, Im1, Im2);
        % postprocess Vectors
        Vectors = PIV_postprocess(PIVParams, Vectors);
        % calculate change in dx's and dy's
        dxChange = max(abs(dxOld - Vectors.dx));
        dyChange = max(abs(dyOld - Vectors.dy));
        maxChange = max(dxChange, dyChange);
        if flag_disp == 1
            display(['maxChange: ' num2str(maxChange)])
        end
        
    else
        %% don't do analysis, but instead load PIV data:
        load(DataFileName)
    end
    
    %% data conversion
    if PIVParams.ConvertData
        Vectors = PIV_conversion(PIVParams, Vectors); 
    end
    
    %% calculate scalar fields
    if PIVParams.CalcScalars
        Vectors = PIV_calcScalarFields(PIVParams, Vectors);
    end
    
    %% Data output
    % generate filename for images and output data file
    if PIVParams.MultipageTiff
        full_filename = [PIVParams.ImageOutputDir PIVParams.files(1).name];
        filename = [full_filename(1:end-4) num2str(i,'%06u')];
    else
        full_filename = [PIVParams.ImageOutputDir PIVParams.files(i).name];
        filename = full_filename(1:end-4);
    end
    
    if PIVParams.SaveData
        % save vectors
        save(DataFileName,'Vectors')
    end
    
    
    %% visualize data
    if PIVParams.DisplayEveryFrame
        PIVParams = PIV_showvectors(PIVParams, Vectors, i, ImOriginal);
        % save image
        if PIVParams.SaveImages
            switch PIVParams.ImageOutputFormat
                case 'png'
                    print(PIVParams.h_fig,'-r300','-dpng', [filename '.png'])
                    % Get rid of the white boarders when there is time
                    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    % ImShow = getframe(PIVParams.h_fig); 
                    % imwrite(ImShow.cdata,[filename '.png'])
                case 'pdf'
                    print(PIVParams.h_fig,'-dpdf', [filename '.pdf'])
            end
        end
    end
    
    %% display global information per frame
    % display maximum displacement
    displacement = sqrt((Vectors.dx).^2 + (Vectors.dy).^2);
    maxDisplacement = max(displacement);
    gridSize = 4 * ceil(maxDisplacement);
    if flag_disp == 1
        disp(['Max displacement: ' num2str(maxDisplacement) ' pixels'])
        disp(['Recommended minimum (initial) grid size: ' num2str(gridSize) ' pixels'])
    end
    
    %% cleanup
    if (i < PIVParams.LastImage)
        clear Vectors
    end
end

toc

%% Did we finish? If not, say why not:
if PIVParams.Abort
    disp(['WARNING: PIV aborted. (' PIVParams.AbortMessage ')'])
    PIV_writeErrorFile(PIVParams)
else
    disp('Finished PIV.')
end
