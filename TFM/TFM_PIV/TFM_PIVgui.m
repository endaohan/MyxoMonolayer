function TFM_PIVgui(varargin)
% PIVgui Brief description of GUI.
%        Comments displayed at the command line in response
%        to the help command.

% (Leave a blank line following the help.)

%%  Initialization tasks

figureHeight = 450;

fh = figure(...
    'MenuBar','none',...
    'Toolbar','none',...
    'Position',[10 10 400 figureHeight],...
    'Color',[0.94 0.94 0.94],...
    'NumberTitle','off',...
    'Name','PIV Parameters',...
    'Visible','off');
movegui(fh,'center')

% defaults:
PIVParams = PIV_setDefaults();

% Save the structure
guidata(fh,PIVParams)

%%  Construct the components

%-----------------------
% Buttons to the different groups of settings
vPos = figureHeight-60;
vSpacing = 40;
buttonWidth = 110;
uicontrol(fh,'Style','pushbutton',...
    'String','File locations...',...
    'Position',[30 vPos buttonWidth 30],...
    'Callback',@FileLocationsCallback);
vPos = vPos - vSpacing;
uicontrol(fh,'Style','pushbutton',...
    'String','Image preprocess...',...
    'Position',[30 vPos buttonWidth 30],...
    'Callback',@ImagePreprocessCallback);
vPos = vPos - vSpacing;
uicontrol(fh,'Style','pushbutton',...
    'String','Grid and ROI...',...
    'Position',[30 vPos buttonWidth 30],...
    'Callback',@GridAndROICallback);
vPos = vPos - vSpacing;
uicontrol(fh,'Style','pushbutton',...
    'String','Correlation...',...
    'Position',[30 vPos buttonWidth 30],...
    'Callback',@CorrelationCallback);
vPos = vPos - vSpacing;
uicontrol(fh,'Style','pushbutton',...
    'String','Validation & filtering...',...
    'Position',[30 vPos buttonWidth 30],...
    'Callback',@ValidationCallback);
vPos = vPos - vSpacing;
uicontrol(fh,'Style','pushbutton',...
    'String','Conversion...',...
    'Position',[30 vPos buttonWidth 30],...
    'Callback',@ConversionCallback);
vPos = vPos - vSpacing;
uicontrol(fh,'Style','pushbutton',...
    'String','Visualization...',...
    'Position',[30 vPos buttonWidth 30],...
    'Callback',@VisualizationCallback)

%----------------------------------
% Tasks to perform on run
vPos = figureHeight-50;
vSpacing = 20;
hPos = 250;
uicontrol(fh,'Style','text',...
    'String','Tasks to perform on run:',...
    'Position',[hPos-20 vPos 130 20]);
vPos = vPos - vSpacing;
DoPIVAnalysis = uicontrol(fh,'Style','checkbox',...
    'String','Do PIV analysis',...
    'Value',PIVParams.DoPIVAnalysis,...
    'Position',[hPos vPos 130 20], ...
    'Callback',@DoPIVAnalysisCallback);
vPos = vPos - vSpacing;
% Convert dx's and dy's to u's and v's in physical units
ConvertData = uicontrol(fh,'Style','checkbox',...
    'String','Convert PIV data',...
    'Value',PIVParams.ConvertData,...
    'Position',[hPos vPos 130 20], ...
    'Callback',@ConvertDataCallback);
vPos = vPos - vSpacing;
% Calculate scalar fields:
CalcScalars = uicontrol(fh,'Style','checkbox',...
    'String','Calculate scalar fields',...
    'Value',PIVParams.CalcScalars,...
    'Position',[hPos vPos 130 20], ...
    'Callback',@CalcScalarsCallback);
vPos = vPos - vSpacing;
% Save data:
SaveData = uicontrol(fh,'Style','checkbox',...
    'String','Save data',...
    'Value',PIVParams.SaveData,...
    'Position',[hPos vPos 130 20], ...
    'Callback',@SaveDataCallback);
vPos = vPos - vSpacing;
% Display the vectors of every frame that is processed or loaded
% (last image is always displayed)
DisplayEveryFrame = uicontrol(fh,'Style','checkbox',...
    'String','Display every frame',...
    'Value',PIVParams.DisplayEveryFrame,...
    'Position',[hPos vPos 130 20], ...
    'Callback',@DisplayEveryFrameCallback);
vPos = vPos - vSpacing;
% Save images (only displayed images are saved)
SaveImages = uicontrol(fh,'Style','checkbox',...
    'String','Save images',...
    'Value',PIVParams.SaveImages,...
    'Position',[hPos vPos 130 20], ...
    'Callback',@SaveImagesCallback);
vPos = vPos - vSpacing - 20;
uicontrol(fh,'Style','pushbutton',...
    'String','Save and run',...
    'Position',[hPos vPos buttonWidth 30],...
    'Callback',@SaveAndRunCallback);


%---------------------------------
% File operations
hPos = 30;
uicontrol(fh,'Style','pushbutton',...
    'String','Save',...
    'Position',[hPos 30 buttonWidth 30],...
    'Callback',@MakeFileCallback);
hPos = hPos + buttonWidth + 10;
uicontrol(fh,'Style','pushbutton',...
    'String','Save as...',...
    'Position',[hPos 30 buttonWidth 30],...
    'Callback',@SaveAsCallback);
hPos = hPos + buttonWidth + 10;
uicontrol(fh,'Style','pushbutton',...
    'String','Open...',...
    'Position',[hPos 30 buttonWidth 30],...
    'Callback',@LoadCallback);





%%  Initialization tasks

set(fh,'Visible','on')

%%  Callbacks
%--------------------------
% Buttons
    function FileLocationsCallback(~, ~)
        PIVgui_FileLocations(fh)
    end
    function ImagePreprocessCallback(~,~)
        PIVgui_ImagePreprocess(fh)
    end
    function GridAndROICallback(~,~)
        PIVgui_GridAndROI(fh)
    end
    function CorrelationCallback(~, ~)
        PIVgui_Correlation(fh)
    end
    function ValidationCallback(~,~)
        PIVgui_Validation(fh)
    end
    function ConversionCallback(~,~)
        PIVgui_Conversion(fh)
    end
    function VisualizationCallback(~,~)
        PIVgui_Visualization(fh)
    end

%-----------------------------
% Tasks to perform
function DoPIVAnalysisCallback(hObject,~)
    PIVParams = guidata(fh);
    PIVParams.DoPIVAnalysis = get(hObject,'Value');
    guidata(fh,PIVParams)
end
function ConvertDataCallback(hObject,~)
    PIVParams = guidata(fh);
    PIVParams.ConvertData = get(hObject,'Value');
    guidata(fh,PIVParams)
end
function CalcScalarsCallback(hObject,~)
    PIVParams = guidata(fh);
    PIVParams.CalcScalars = get(hObject,'Value');
    guidata(fh,PIVParams)
end
function SaveDataCallback(hObject,~)
    PIVParams = guidata(fh);
    PIVParams.SaveData = get(hObject,'Value');
    guidata(fh,PIVParams)
end
function DisplayEveryFrameCallback(hObject,~)
    PIVParams = guidata(fh);
    PIVParams.DisplayEveryFrame = get(hObject,'Value');
    if ~PIVParams.DisplayEveryFrame
        PIVParams.SaveImages = false;
        set(SaveImages,'Value',false)
    end
    guidata(fh,PIVParams)
end
function SaveImagesCallback(hObject,~)
    PIVParams = guidata(fh);
    PIVParams.SaveImages = get(hObject,'Value');
    if PIVParams.SaveImages
        PIVParams.DisplayEveryFrame = true;
        set(DisplayEveryFrame,'Value',true)
    end
    guidata(fh,PIVParams)
end
function SaveAndRunCallback(~, ~)
    if isempty(PIVParams.CurrentParametersFile)
        saveAs
    else
        makeFile(PIVParams.CurrentParametersFile)
    end
    if ~isempty(PIVParams.CurrentParametersFile)
        TFM_PIV(PIVParams.CurrentParametersFile)
    end
end


%--------------------------------------
% Parameters file saving, loading, recovering
    function MakeFileCallback(~,~)
        if isempty(PIVParams.CurrentParametersFile)
            saveAs
        else
            makeFile(PIVParams.CurrentParametersFile)
        end
    end
    function SaveAsCallback(~,~)
        saveAs
    end
    function LoadCallback(~,~)
        % let user select a file
        [FileName,PathName,~] = uigetfile([PIVParams.CurrentDir '/*.txt']);
        % check if file selection was succesfull
        if ischar(FileName)
            % read the selected file
            PIVParams = readPIVParameters([PathName FileName]);
            % set the current directory to the selected path
            PIVParams.CurrentDir = PathName;
            % empty CurrentParametersFile to prevent accidentally
            % overwriting
            PIVParams.CurrentParametersFile = '';
            % set checkboxes in main window to the correct values
            set(DoPIVAnalysis,'Value',PIVParams.DoPIVAnalysis)
            set(ConvertData,'Value',PIVParams.ConvertData)
            set(CalcScalars,'Value',PIVParams.CalcScalars)
            set(SaveData,'Value',PIVParams.SaveData)
            set(DisplayEveryFrame,'Value',PIVParams.DisplayEveryFrame)
            set(SaveImages,'Value',PIVParams.SaveImages)
            % save the structure
            guidata(fh,PIVParams)
        end
    end

%%  Utility functions
    function makeFile(fileName)
        % format: each line consists of a string (%s) and newline command (\r\n)
        format = '%s\r\n';
        % make a true / false cell for easy entry of logical values
        logicalEntry = {'false', 'true'};
        % load the data
        PIVParams = guidata(fh);
        
        % open output file
        fid = fopen(fileName,'w');
        % write header to file
        fprintf(fid, format, '#PIV-Parameters');
        fprintf(fid, '\r\n');
        fprintf(fid, format, '% PIV Parameters file generated by PIVgui');
        fprintf(fid, '\r\n');
        fprintf(fid, format, ['% generated on: ' datestr(now,'mmmm dd, yyyy HH:MM:SS AM')]);
        fprintf(fid, '\r\n');
        fprintf(fid, '\r\n');
        fprintf(fid, '\r\n');
        
        %% -----------------------------
        fprintf(fid, format, '%%%%%%%%%%%%%');
        fprintf(fid, format, '%%% Tasks %%%');
        fprintf(fid, format, '%%%%%%%%%%%%%');
        fprintf(fid, format, ['#DoPIVAnalysis: ' char(logicalEntry(PIVParams.DoPIVAnalysis+1))]);
        fprintf(fid, format, ['#ConvertData: ' char(logicalEntry(PIVParams.ConvertData+1))]);
        fprintf(fid, format, ['#CalcScalars: ' char(logicalEntry(PIVParams.CalcScalars+1))]);
        fprintf(fid, format, ['#SaveData: ' char(logicalEntry(PIVParams.SaveData+1))]);
        fprintf(fid, format, ['#DisplayEveryFrame: ' char(logicalEntry(PIVParams.DisplayEveryFrame+1))]);
        fprintf(fid, format, ['#SaveImages: ' char(logicalEntry(PIVParams.SaveImages+1))]);
        fprintf(fid, '\r\n');
        
        %% ------------------------------
        % File locations and names
        fprintf(fid, format, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fid, format, '%%% File locations and names %%%');
        fprintf(fid, format, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fid, format, ['#ImageFileName: ' PIVParams.ImageFileName]);
        fprintf(fid, format, ['#FileExtension: ' PIVParams.FileExtension]);
        fprintf(fid, format, ['#MultipageTiff: ' char(logicalEntry(PIVParams.MultipageTiff+1))]);
        fprintf(fid, format, ['#ImageBitDepth: ' num2str(PIVParams.ImageBitDepth)]);
        fprintf(fid, format, ['#Directory: ' PIVParams.Directory]);
        fprintf(fid, format, ['#ImageOutputDir: ' PIVParams.ImageOutputDir]);
        fprintf(fid, format, ['#ImageOutputFormat: ' PIVParams.ImageOutputFormat]);
        fprintf(fid, format, ['#DataDir: ' PIVParams.DataDir]);
        fprintf(fid, format, ['#FirstImage: ' num2str(PIVParams.FirstImage)]);
        fprintf(fid, format, ['#LastImage: ' num2str(PIVParams.LastImage)]);
        fprintf(fid, format, ['#SkipImages: ' num2str(PIVParams.SkipImages)]);
        fprintf(fid, '\r\n');
        
        %% ------------------------------
        % Image pre-processing
        fprintf(fid, format, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fid, format, '%%% Image pre-processing %%%');
        fprintf(fid, format, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fid, format, ['#UseMinMax: ' char(logicalEntry(PIVParams.UseMinMax+1))]);
        fprintf(fid, format, ['#MinMaxKernel: ' num2str(PIVParams.MinMaxKernel)]);
        fprintf(fid, format, ['#UseUnsharpMask: ' char(logicalEntry(PIVParams.UseUnsharpMask+1))]);
        fprintf(fid, format, ['#UnsharpMaskAmount: ' num2str(PIVParams.UnsharpMaskAmount)]);
        fprintf(fid, format, ['#UseGaussianBlur: ' char(logicalEntry(PIVParams.UseGaussianBlur+1))]);
        fprintf(fid, format, ['#GaussianBlurRadius: ' num2str(PIVParams.GaussianBlurRadius)]);
        fprintf(fid, '\r\n');
        
        %% ------------------------------
        % Grid and ROI
        fprintf(fid, format, '%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fid, format, '%%% Grid size and ROI %%%');
        fprintf(fid, format, '%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fid, format, ['#RotationAngle: ' num2str(PIVParams.RotationAngle)]);
        fprintf(fid, format, ['#GridSize: ' num2str(PIVParams.GridSize)]);
        fprintf(fid, format, ['#Overlap: ' num2str(PIVParams.Overlap)]);
        fprintf(fid, format, ['#UseROI: ' char(logicalEntry(PIVParams.UseROI+1))]);
        fprintf(fid, format, ['#CropImageToROI: ' char(logicalEntry(PIVParams.CropImageToROI+1))]);
        fprintf(fid, format, ['#ROILeftEdge: ' num2str(PIVParams.ROILeftEdge)]);
        fprintf(fid, format, ['#ROIRightEdge: ' num2str(PIVParams.ROIRightEdge)]);
        fprintf(fid, format, ['#ROITopEdge: ' num2str(PIVParams.ROITopEdge)]);
        fprintf(fid, format, ['#ROIBottomEdge: ' num2str(PIVParams.ROIBottomEdge)]);
        fprintf(fid, format, ['#UseMovingROI: ' char(logicalEntry(PIVParams.UseMovingROI+1))]);
        fprintf(fid, format, ['#MovingROILeftStart: ' num2str(PIVParams.MovingROILeftStart)]);
        fprintf(fid, format, ['#MovingROIRightStart: ' num2str(PIVParams.MovingROIRightStart)]);
        fprintf(fid, format, ['#MovingROIHorizontalVelocity: ' num2str(PIVParams.MovingROIHorizontalVelocity)]);
        fprintf(fid, format, ['#MovingROITopStart: ' num2str(PIVParams.MovingROITopStart)]);
        fprintf(fid, format, ['#MovingROIBottomStart: ' num2str(PIVParams.MovingROIBottomStart)]);
        fprintf(fid, format, ['#MovingROIVerticalVelocity: ' num2str(PIVParams.MovingROIVerticalVelocity)]);
        fprintf(fid, '\r\n');
        
        %% ------------------------------
        % Correlation
        fprintf(fid, format, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fid, format, '%%% Correlation settings %%%');
        fprintf(fid, format, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fid, format, ['#Method: ' PIVParams.Method]);
        fprintf(fid, format, ['#NumPasses: ' num2str(PIVParams.NumPasses)]);
        fprintf(fid, format, ['#UseGridRefinement: ' char(logicalEntry(PIVParams.UseGridRefinement+1))]);
        fprintf(fid, format, ['#InitialGridSize: ' num2str(PIVParams.InitialGridSize)]);
        fprintf(fid, format, ['#CorrelationStep: ' num2str(PIVParams.CorrelationStep)]);
        fprintf(fid, '\r\n');
        
        %% ------------------------------
        % Validation and filtering
        fprintf(fid, format, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fid, format, '%%% Validation and filtering %%%');
        fprintf(fid, format, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fid, format, ['#numValidationPasses: ' num2str(PIVParams.numValidationPasses)]);
        fprintf(fid, format, ['#UseMaxDisplacement: ' char(logicalEntry(PIVParams.UseMaxDisplacement+1))]);
        fprintf(fid, format, ['#MaxDisplacement: ' num2str(PIVParams.MaxDisplacement)]);
        fprintf(fid, format, ['#UseMinCorrelationCoefficient: ' char(logicalEntry(PIVParams.UseMinCorrelationCoefficient+1))]);
        fprintf(fid, format, ['#MinCorrelationCoefficient: ' num2str(PIVParams.MinCorrelationCoefficient)]);
        fprintf(fid, format, ['#UseNormalizedMedianTest: ' char(logicalEntry(PIVParams.UseNormalizedMedianTest+1))]);
        fprintf(fid, format, ['#MedianThreshold: ' num2str(PIVParams.MedianThreshold)]);
        fprintf(fid, format, ['#UseSecondPeak: ' char(logicalEntry(PIVParams.UseSecondPeak+1))]);
        fprintf(fid, format, ['#UseInterpolation: ' char(logicalEntry(PIVParams.UseInterpolation+1))]);
        fprintf(fid, format, ['#UseDataSmoothing: ' char(logicalEntry(PIVParams.UseDataSmoothing+1))]);
        fprintf(fid, format, ['#SmoothingFilterSize: ' num2str(PIVParams.SmoothingFilterSize)]);
        fprintf(fid, format, ['#SmoothingFilterType: ' PIVParams.SmoothingFilterType]);
        fprintf(fid, '\r\n');
        
        %% ------------------------------
        % Conversion
        fprintf(fid, format, '%%%%%%%%%%%%%%%%%%');
        fprintf(fid, format, '%%% Conversion %%%');
        fprintf(fid, format, '%%%%%%%%%%%%%%%%%%');
        fprintf(fid, format, ['#FPS: ' num2str(PIVParams.FPS)]);
        fprintf(fid, format, ['#Scaling: ' num2str(PIVParams.Scaling)]);
        fprintf(fid, format, ['#TimeUnit: ' PIVParams.TimeUnit]);
        fprintf(fid, format, ['#LengthUnit: ' PIVParams.LengthUnit]);
        fprintf(fid, format, ['#VelocityUnit: ' PIVParams.VelocityUnit]);
        fprintf(fid, '\r\n');
        
        %% ------------------------------
        % Visualization
        fprintf(fid, format, '%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fid, format, '%%% Visualization %%%');
        fprintf(fid, format, '%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fid, format, ['#ShowVectors: ' char(logicalEntry(PIVParams.ShowVectors+1))]);
        fprintf(fid, format, ['#ColorCodeVectors: ' char(logicalEntry(PIVParams.ColorCodeVectors+1))]);
        fprintf(fid, format, ['#VectorColorValidation: ' char(logicalEntry(PIVParams.VectorColorValidation+1))]);
        fprintf(fid, format, ['#VectorColormap: ' PIVParams.VectorColormap]);
        if ischar(PIVParams.VectorColor)
            fprintf(fid, format, ['#VectorColor: ' PIVParams.VectorColor]);
        else
            fprintf(fid, format, ['#VectorColor: [' num2str(PIVParams.VectorColor) ']']);
        end
        fprintf(fid, format, ['#VectorScale: ' num2str(PIVParams.VectorScale)]);
        fprintf(fid, format, ['#ShowImage: ' char(logicalEntry(PIVParams.ShowImage+1))]);
        fprintf(fid, format, ['#ShowGrid: ' char(logicalEntry(PIVParams.ShowGrid+1))]);
        fprintf(fid, format, ['#ShowContour: ' char(logicalEntry(PIVParams.ShowContour+1))]);
        fprintf(fid, format, ['#ContourScalar: ' PIVParams.ContourScalar]);
        fprintf(fid, format, ['#ContourColormap: ' PIVParams.ContourColormap]);
        fprintf(fid, format, ['#ContourSteps: ' num2str(PIVParams.ContourSteps)]);
        fprintf(fid, format, ['#ScalarMinValue: ' num2str(PIVParams.ScalarMinValue)]);
        fprintf(fid, format, ['#ScalarMaxValue: ' num2str(PIVParams.ScalarMaxValue)]);
        fprintf(fid, format, ['#ScalarAutoScale: ' char(logicalEntry(PIVParams.ScalarAutoScale+1))]);
        fprintf(fid, format, ['#ScalarPrefactor: ' num2str(PIVParams.ScalarPrefactor)]);
        fprintf(fid, format, ['#DisplayScalarInfo: ' char(logicalEntry(PIVParams.DisplayScalarInfo+1))]);
        fprintf(fid, format, ['#OutputImageWidth: ' num2str(PIVParams.OutputImageWidth)]);
        fprintf(fid, format, ['#OutputImageHeight: ' num2str(PIVParams.OutputImageHeight)]);
        fprintf(fid, format, ['#FontName: ' PIVParams.FontName]);
        fprintf(fid, format, ['#FontSize: ' num2str(PIVParams.FontSize)]);
        fprintf(fid, format, ['#ScalarLabel: ' PIVParams.ScalarLabel]);
        
        % close file
        fclose(fid);
    end
    function saveAs
        % get a file name and location
        [FileName,PathName] = uiputfile([PIVParams.CurrentDir '/parameters.txt']);
        if ischar(FileName)
            % load the structure
            PIVParams = guidata(fh);
            % adjust the data
            PIVParams.CurrentParametersFile = [PathName FileName];
            PIVParams.CurrentDir = PathName;
            % save the structure
            guidata(fh,PIVParams)
            % save data to the selected file
            makeFile(PIVParams.CurrentParametersFile)
        end
    end
end