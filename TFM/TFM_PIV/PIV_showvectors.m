% Modified by Endao Han, 2019/9/3

function PIVParams = PIV_showvectors(PIVParams, Vectors, imageNumber, Image)

%
% PIV_showvectors - Display the image, contour and/or velocity vectors
% (This function is called by the function PIV)
%

% %% close the old figure (if it exists) and create a new figure window
% if PIVParams.h_fig
%     close(PIVParams.h_fig)
% end
% PIVParams.h_fig = figure;

%% clear the old figure (if it exists) and create a new figure window
if (PIVParams.h_fig ~=0)
    figure(PIVParams.h_fig)
    clf
else
    width = PIVParams.OutputImageWidth;
    height = PIVParams.OutputImageHeight;
    PIVParams.h_fig = figure;
    set(gcf,'PaperSize',[width height])
    set(gcf,'PaperPosition',[0 0 width height])
%     set(gcf,'PaperSize',[8 4.5])
%     set(gcf,'PaperPosition',[0 0 8 4.5])
end

%% some other settings
fontName = PIVParams.FontName;
fontSize = PIVParams.FontSize;

%% read the image if needed
if nargin<4
    if PIVParams.MultipageTiff
        % read image from multipage tiff file
        Im1 = imread([PIVParams.Directory PIVParams.files(1).name], imageNumber);
    else
        % read image
        Im1 = imread([PIVParams.Directory PIVParams.files(imageNumber).name]);
    end
    % preprocess images
    [Image, ~] = PIV_prepocess(PIVParams, Vectors, Im1, Im1);
end

%% create axes
% imageAxes = axes;
axis equal
% set limits
[ImageHeight, ImageWidth] = size(Image);
xlim([1 ImageWidth])
ylim([1 ImageHeight])
hold on
set(gca,'YDir','reverse', ...
    'XTick',[], ...
    'YTick',[])
box on

%% lowest layer: the image
% display the image:
if PIVParams.ShowImage
    % gamma correction
    % gamma = 0.5; 
    gamma = 1; 
    
    % Added by Endao, to use rgb images
    % I = double(rgb2gray(Image)) / 255;
    
    I = double(Image) / 255;
    I = I.^gamma;
    I = uint8(I * 255);
    subimage(I,gray(256))
%     subimage(uint8(0.25*Image),gray(64))
%     subimage(uint8(Image),gray(256))
    %,'InitialMagnification','fit')
%     imshow(uint8(Image),'InitialMagnification','fit')
    hold on
%     colormap(imageAxes,gray)
end

%% next layer: contour of scalar value
if PIVParams.ShowContour
    % default units for the colorbar:
    units = '';
    % Get the coordinates
    y = Vectors.y;
    x = Vectors.x;
    % select a scalar to show:
    switch PIVParams.ContourScalar
        case 'magnitude'
            z = Vectors.magnitude;
            units = PIVParams.VelocityUnit;
        case 'signalToNoiseRatio'
            z = Vectors.signalToNoiseRatio;
        case 'correlationCoefficient'
            z = Vectors.correlationCoefficient;
        case 'vorticity'
            z = Vectors.vorticity;
            units = [PIVParams.TimeUnit '^{-1}'];
        case 'shearStrain'
            z = Vectors.shearStrain;
            units = [PIVParams.TimeUnit '^{-2}'];
        case 'normalStrain'
            z = Vectors.normalStrain;
            units = ['\itdu\rm/\itdx\rm + \itdv\rm/\itdy\rm (' PIVParams.TimeUnit '^{-1})'];
        case 'dudx'
            z = Vectors.dudx;
            units = ['du/dx (' PIVParams.TimeUnit '^{-1})'];
        case 'dudy'
            z = Vectors.dudy;
            units = [PIVParams.TimeUnit '^{-1}'];
        case 'dissipation'
            z = (Vectors.dudy.^2 + Vectors.dvdx.^2);
        case 'dvdx'
            z = Vectors.dvdx;
            units = ['dv/dx (' PIVParams.TimeUnit '^{-1})'];
        case 'dvdy'
            z = Vectors.dvdy;
            units = ['\itdu_y\rm/\itdy\rm (' PIVParams.TimeUnit '^{-1})'];
        case 'xComponent'
            z = Vectors.u;
            units = ['\itu_x\rm (' (PIVParams.VelocityUnit) ')'];
        case 'yComponent'
            z = Vectors.v;
            units = PIVParams.VelocityUnit;
            
    end
    
    % apply custom units
    if ~isempty(PIVParams.ScalarLabel)
        units = PIVParams.ScalarLabel;
    end
    
    % apply prefactor for scalar field
    zPrefactor = PIVParams.ScalarPrefactor;
    z = zPrefactor * z;
    % take absolute value
    %z = abs(z);
    % get dimensions of vector field
    width = find(y==y(1),1,'last');
%     width = i-1;
    height = length(x) / width;
    %reshape the vectors into matrices:
    x = reshape(x,width,height)';
    y = reshape(y,width,height)';
    z = reshape(z,width,height)';
    
    if PIVParams.DisplayScalarInfo
        display(['Contour:       ' PIVParams.ContourScalar])
        display(['Minimum value: ' num2str(min(min(z)))])
        display(['Maximum value: ' num2str(max(max(z)))])
        display(['Mean value:    ' num2str(mean(mean(z)))])
    end
    
    % rescale scalar field
    if PIVParams.ScalarAutoScale
        minZ = min(min(z));
        maxZ = max(max(z));
    else
        minZ = PIVParams.ScalarMinValue;
        maxZ = PIVParams.ScalarMaxValue;
    end
    
    % draw lines
    if PIVParams.ContourLines
        LineStyle = '-';
    else
        LineStyle = 'none';
    end
    
    % plot the contour
    LevelList = linspace(minZ,maxZ,PIVParams.ContourSteps+1);
    LevelList = [min(z(:)) LevelList(2:end-1) max(z(:))];
    z(z==0)=NaN;
%     figure
%     contourf(x,y,z,'LevelList',LevelList,'LineStyle',LineStyle);
%     pause
%     close
    contourf(x,y,z,'LevelList',LevelList,'LineStyle',LineStyle);
    %     contourf(x,y,z,PIVParams.ContourSteps-1,'LineStyle',LineStyle);
    
    % make nice colors
    cmap = [];
    switch PIVParams.ContourColormap
        case 'BlueWhiteRed'
            % this is just to be compatible with old parameter files
            cmap = 'PIVBlueWhiteRed';
            eval(['colormap(gca,' cmap '(' num2str(PIVParams.ContourSteps) '))'])
        case 'GreenYellowRed'
            load PIVColors_GreenYellowRed.mat
            colormap(gca,cmap)
        otherwise
            cmap = PIVParams.ContourColormap;
            eval(['colormap(gca,' cmap '(' num2str(PIVParams.ContourSteps) '))'])
    end
    
    
    % adjust the scale of the colors
    set(gca,'CLim',[minZ maxZ])
    h = colorbar;
    xlabel(h,units,'FontName',fontName);
    set(gca,'FontName',fontName)
    % xlabel(h,units,'FontSize',20);
    % set(h,'FontSize',20)
    
end


%% next layer: grid
if PIVParams.ShowGrid
%     % display the grid
%     [ImageHeight, ImageWidth] = size(Image);
%     % vertical grid lines
%     i_max = floor(ImageWidth/PIVParams.GridSize);
%     for i=1:i_max
%         line([1,1]*i*PIVParams.GridSize,[1, ImageHeight])
%     end
%     % horizontal grid lines
%     i_max = floor(ImageHeight/PIVParams.GridSize);
%     for i=1:i_max
%         line([1, ImageWidth], [1,1]*i*PIVParams.GridSize)
%     end
    
    % get the coordinates
    x = Vectors.x;
    y = Vectors.y;
    grid = PIVParams.GridSize;
    minX = round(x(1)-grid/2);
    maxX = round(max(x)+grid/2);
    minY = round(y(1)-grid/2);
    maxY = round(max(y)+grid/2);
    for i=minX:grid:maxX
        line([i,i],[minY,maxY])
    end
    for i=minY:grid:maxY
        line([minX,maxX],[i,i])
    end
    
    
end

%% top layer: the vectors
if PIVParams.ShowVectors
    % display settings
    color = PIVParams.VectorColor;
    % get the coordinates
    x = Vectors.x;
    y = Vectors.y;
    % the velocities
    u = Vectors.u;
    v = Vectors.v;
    
    
    scale = 0; % this disables the automatic MATLAB scaling (because we don't want MATLAB to
    % autoscale our vectors.
    
    if PIVParams.VectorColorValidation
        %%% Color code the vectors according to the validation flag
        colorValid = 'green';
        colorInvalid = 'red';
        colorReplaced = [1.0 0.5 0.5];
        
        LineWidth = 0.5;
        
        u = u * PIVParams.VectorScale;
        v = v * PIVParams.VectorScale;
        
        x(v==0)=NaN;
        y(u==0)=NaN;

        
        % valid vectors:
        i = Vectors.validationFlag==0;
        quiver(x(i), y(i), u(i), -v(i), scale, 'color', colorValid, 'LineWidth', LineWidth)
        
        % invalid vectors:
        i = Vectors.validationFlag==1;
        quiver(x(i), y(i), u(i), -v(i), scale, 'color', colorInvalid, 'LineWidth', LineWidth)
        
        % replaced vectors:
        i = Vectors.validationFlag==2;
        quiver(x(i), y(i), u(i), -v(i), scale, 'color', colorReplaced, 'LineWidth', LineWidth)
        
    else
        
        % find vectors that we shouldn't plot:
        validIndices = logical(Vectors.validationFlag~=1);
        % do not plot invalid vectors:
        x = x(validIndices);
        y = y(validIndices);
        u = u(validIndices);
        v = v(validIndices);
        
        
        
        
        % do the actual displaying of the vectors:
        if PIVParams.ColorCodeVectors
            % create a nice color coded vector plot
            PIV_colorQuiver(PIVParams, Vectors);
        else
            u = u * PIVParams.VectorScale;
            v = v * PIVParams.VectorScale;
            quiver(x, y, u, -v, scale, 'color', color)
        end
    end
end

%% render the figure
drawnow