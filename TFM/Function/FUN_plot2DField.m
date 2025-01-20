% *********************************************************************** %
% FUN_plot2DField( U, V, x, y, imgPolar, imgDir, cMap, n_box, vecRatio,flag_bar)
% 
% Function to plot 2D vector field. 
% Colormap represents the magnitude of the field, and arrows show the
% direction and magnitude of the vectors. 
% Can also coarse-grain vectors to make them more visible. 
% Written by Endao Han, V4, 2021/1/17. 
% *********************************************************************** %

function FUN_plot2DField( U, V, x, y, imgPolar, imgDir, cMap, n_box, vecRatio,flag_bar )

% set(hd,'units','inch','position',[1,0,10,5])

if numel(U) ~= numel(V)
    disp('ERROR: FUN_plot2DField: input matrices do not have the same size. ')
    return
end

% Magnitude of the field
mag = sqrt( U.^2+V.^2 ); 
n = size(mag); 
cmax = max( mag,[],'all' ); 

% Polarity of the image: 1 - all positive; 0 - on both sides
if isempty(imgPolar)
    imgPolar = 0; 
end

switch imgPolar
    case 0
        % Set the boundaries of colormap
        cRange = [-cmax cmax]; 
        % Default colormap is BlueWhiteRed
        if isempty(cMap)
            cMap = 'BlueWhiteRed'; 
        end
    case 1
        % Set the boundaries of colormap
        cRange = [0 cmax]; 
        % Default colormap is BlueWhiteRed
        if isempty(cMap)
            cMap = 'parula'; 
        end
end

% Default image direction
if isempty(imgDir)
    % imgDir = 'normal'; 
    imgDir = 'reverse'; 
end

% Default vector ratio
if ~exist('vecRatio','var')
    vecRatio = 20/cmax; 
end

% Default flag for color bar
if nargin == 9
    flag_bar = 1; 
end
    


% Set default x, y positions
if ( isempty(x) || isempty(y) )
    x = linspace( 1,n(2),n(2) ); 
    y = linspace( 1,n(1),n(1) ); 
end

% Coarse-grain vectors
if n_box ~= 0 
    N_box = floor(n/n_box); 

    U_crop = U(1:N_box(1)*n_box,1:N_box(2)*n_box); 
    V_crop = V(1:N_box(1)*n_box,1:N_box(2)*n_box); 

    U_box = imresize( U_crop,N_box,'box' ); 
    V_box = imresize( V_crop,N_box,'box' ); 
    
    i_box = round( ( linspace(1,N_box(2),N_box(2))-0.5 )*n_box ); 
    j_box = round( ( linspace(1,N_box(1),N_box(1))-0.5 )*n_box ); 
     
    x_box = x(i_box); 
    y_box = y(j_box); 
end


% Make plot
imagesc(x,y,mag)
hold on
if n_box > 1
    if strcmp(imgDir,'reverse')
        quiver( x_box,y_box, U_box*vecRatio,-V_box*vecRatio, ... 
         'k',  'AutoScale','off'); 
    elseif strcmp(imgDir,'normal')
         quiver( x_box,y_box, U_box*vecRatio,V_box*vecRatio, ... 
             'k',  'AutoScale','off'); 
    else
        disp('EH ERROR: Wrong variable imgDir!!')
    end
else
    if strcmp(imgDir,'reverse')
        quiver( x,y, U*vecRatio,-V*vecRatio, ... 
             'k',  'AutoScale','off'); 
    elseif strcmp(imgDir,'normal')
        quiver( x,y, U*vecRatio,V*vecRatio, ... 
             'k',  'AutoScale','off'); 
    else
        disp('EH ERROR: Wrong variable imgDir!!')
    end      
end
hold off
ax = gca;
ax.YDir = imgDir; 
axis equal
axis( [min(x) max(x) min(y) max(y)] ); 

colormap(cMap); 
caxis(cRange); 
if flag_bar == 1
    colorbar; 
end


end



