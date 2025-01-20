function PIVgui_GridAndROI(main_handle)
% PIVgui Brief description of GUI.
%        Comments displayed at the command line in response
%        to the help command.

% (Leave a blank line following the help.)

%%  Initialization tasks

figureHeight = 500;
figureWidth = 600;

fh = figure(...
    'MenuBar','none',...
    'Toolbar','none',...
    'Position',[10 10 figureWidth figureHeight],...
    'Color',[0.94 0.94 0.94],...
    'NumberTitle','off',...
    'Name','Grid and ROI',...
    'Visible','off',...
    'WindowStyle','modal');
movegui(fh,'center')
set(fh,'CloseRequestFcn',@CancelCallback)

% load data
PIVParams = guidata(main_handle);

%%  Construct the components
vPos = figureHeight-50;
vSpacing = 40;
hPos = 20;
hSpacing = 170;
textWidth = hSpacing-10;
editWidth = 80;
controlHeight = 25;

uicontrol(fh,'Style','text',...
    'String','Grid size',...
    'HorizontalAlignment','Left',...
    'Position',[hPos vPos-5 textWidth controlHeight]);
GridSize = uicontrol(fh,'Style','edit',...
    'String',num2str(PIVParams.GridSize),...
    'Position',[hPos+hSpacing vPos editWidth controlHeight]);
vPos = vPos - vSpacing;
uicontrol(fh,'Style','text',...
    'String','Interrogation window overlap',...
    'HorizontalAlignment','Left',...
    'Position',[hPos vPos-5 textWidth controlHeight]);
Overlap = uicontrol(fh,'Style','edit',...
    'String',num2str(PIVParams.Overlap),...
    'Position',[hPos+hSpacing vPos editWidth controlHeight]);


vPos = vPos - vSpacing - 50;
UseROI = uicontrol(fh,'Style','checkbox',...
    'String','Use ROI',...
    'Value',PIVParams.UseROI,...
    'Position',[hPos+10 vPos 200 controlHeight]);
vPos = vPos - vSpacing;
CropImageToROI = uicontrol(fh,'Style','checkbox',...
    'String','Crop image to ROI',...
    'Value',PIVParams.CropImageToROI,...
    'Position',[hPos+10 vPos 200 controlHeight]);
vPos = vPos + vSpacing / 2;
hPos = hPos + 350;
uicontrol(fh,'Style','pushbutton',...
    'String','Fit ROI to image',...
    'Position',[hPos vPos-5 110 30],...
    'Callback',@FitROIToImageCallback);
hPos = hPos - 130;
editWidth = 50;
uicontrol(fh,'Style','text',...
    'String','ROI',...
    'HorizontalAlignment','Center',...
    'Position',[hPos vPos-5 editWidth controlHeight]);
hPos = hPos - editWidth;
ROILeftEdge = uicontrol(fh,'Style','edit',...
    'String',num2str(PIVParams.ROILeftEdge),...
    'Position',[hPos-5 vPos editWidth controlHeight]);
hPos = hPos + editWidth + editWidth;
ROIRightEdge = uicontrol(fh,'Style','edit',...
    'String',num2str(PIVParams.ROIRightEdge),...
    'Position',[hPos+5 vPos editWidth controlHeight]);
hPos = hPos - editWidth;
vPos = vPos + controlHeight;
ROITopEdge = uicontrol(fh,'Style','edit',...
    'String',num2str(PIVParams.ROITopEdge),...
    'Position',[hPos vPos+5 editWidth controlHeight]);
vPos = vPos - 2*controlHeight;
ROIBottomEdge = uicontrol(fh,'Style','edit',...
    'String',num2str(PIVParams.ROIBottomEdge),...
    'Position',[hPos vPos-5 editWidth controlHeight]);




vPos = vPos - 100;
hPos = 20;
hSpacing = 200;
textWidth = hSpacing-10;
UseMovingROI = uicontrol(fh,'Style','checkbox',...
    'String','Use moving ROI',...
    'Value',PIVParams.UseMovingROI,...
    'Position',[hPos+10 vPos 200 controlHeight]);
vPos = vPos - vSpacing;
uicontrol(fh,'Style','text',...
    'String','Horizontal velocity (pixels/frams):',...
    'HorizontalAlignment','Left',...
    'Position',[hPos vPos-5 textWidth controlHeight]);
MovingROIHorizontalVelocity = uicontrol(fh,'Style','edit',...
    'String',num2str(PIVParams.MovingROIHorizontalVelocity),...
    'Position',[hPos+hSpacing vPos editWidth controlHeight]);
vPos = vPos - vSpacing;
uicontrol(fh,'Style','text',...
    'String','Vertical velocity (pixels/frams):',...
    'HorizontalAlignment','Left',...
    'Position',[hPos vPos-5 textWidth controlHeight]);
MovingROIVerticalVelocity = uicontrol(fh,'Style','edit',...
    'String',num2str(PIVParams.MovingROIVerticalVelocity),...
    'Position',[hPos+hSpacing vPos editWidth controlHeight]);







%%% OK and cancel
uicontrol(fh,'Style','pushbutton',...
    'String','OK',...
    'Position',[30 30 110 30],...
    'Callback',@ApplyCallback);
uicontrol(fh,'Style','pushbutton',...
    'String','Cancel',...
    'Position',[170 30 110 30],...
    'Callback',@CancelCallback);

%%  Initialization tasks

set(fh,'Visible','on')

%%  Callbacks
    function FitROIToImageCallback(~,~)
        % path for image files:
        path = [PIVParams.Directory PIVParams.ImageFileName '*.' PIVParams.FileExtension];
        % generate file list:
        files = dir(path);
        % get the filename (including full path) of the first image:
        filename = [PIVParams.Directory files(1).name];
        % get info from this file
        info = imfinfo(filename);
        % use this info to fill ROI
        set(ROILeftEdge,'String','1');
        set(ROIRightEdge,'String',num2str(info(1).Width));
        set(ROITopEdge,'String','1');
        set(ROIBottomEdge,'String',num2str(info(1).Height));
    end
    function CancelCallback(~,~)
        delete(fh)
    end
    function ApplyCallback(~,~)
        % load the structure
        PIVParams = guidata(main_handle);

        %%% adjust the data
        PIVParams.GridSize = uint16(str2double(get(GridSize,'String')));
        PIVParams.Overlap = uint16(str2double(get(Overlap,'String')));
        PIVParams.UseROI = get(UseROI,'Value');
        PIVParams.CropImageToROI = get(CropImageToROI,'Value');
        PIVParams.ROILeftEdge = uint16(str2double(get(ROILeftEdge,'String')));
        PIVParams.ROIRightEdge = uint16(str2double(get(ROIRightEdge,'String')));
        PIVParams.ROITopEdge = uint16(str2double(get(ROITopEdge,'String')));
        PIVParams.ROIBottomEdge = uint16(str2double(get(ROIBottomEdge,'String')));
        if get(UseMovingROI,'Value')
            PIVParams.UseMovingROI = true;
            PIVParams.MovingROILeftStart = PIVParams.ROILeftEdge;
            PIVParams.MovingROIRightStart = PIVParams.ROIRightEdge;
            PIVParams.MovingROITopStart = PIVParams.ROITopEdge;
            PIVParams.MovingROIBottomStart = PIVParams.ROIBottomEdge;
        end
        PIVParams.MovingROIHorizontalVelocity = str2double(get(MovingROIHorizontalVelocity,'String'));
        PIVParams.MovingROIVerticalVelocity = str2double(get(MovingROIVerticalVelocity,'String'));
        
        %%% save the structure
        guidata(main_handle,PIVParams)
        delete(fh)
    end

%%  Utility functions

end
