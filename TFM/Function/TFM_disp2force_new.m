%*************************************************************************
% TFM_disp2force( folder, fileName, Parameter, chl )
% 
% Based on TFM_disp2force
% Get the force map from the displacement field
% Written by Endao Han, V8, 2021/2/9
% 
% Make this a function that can be called by TFM_MainFUN.m
% All the old versions of this code used the wrong Green's function
% generator, so stop using the versions before V7. 
% When it is finished, save it as FUN_Disp2Force in the TFM_analysis folder. 
% In this version, overlay director field. 
%*************************************************************************

function Trac = TFM_disp2force_new( folder, fileName, Parameter, chl )

% For debugging
%{
clear 
close all

% folder = 'D:\DATA_Confocal\20210126'; 
% fileName = 'PAA430Pa_Beads100nm_noCover_Coat_TimeSeries_Depth_60x_1x_1'; 

folder = 'D:\DATA_Confocal\20210205'; 
fileName = 'PAA430Pa_Beads100nm_noCover_Coat_TimeSeries_Depth_60x_1x_4'; 


% Load parameters
load([folder '\' fileName '\Parameters.mat']); 

% Channel to analyze
chl = 1; 

% Correction coefficient for reconstruction
Parameter.lambda = 1/1E9; 
% lambda = 0; 
%}

%% Load Disp data
route = [folder '\' fileName]; 
routeIn = [route '\Beads_' num2str(chl) '\DispField']; 
load([routeIn '\Disp.mat'],'Disp');


%% Set parameters
% ======================================================================= %
% Location of the folder with the Green's functions
folderG = pwd; 

% Grid size
% d = 10; 
% delta = 0.1; 

% For Myxo
qScale = 0.05; 
% For Pseudomonas
% qScale = 0.004; 
% Make video? Yes - 1; No -2. 
flag_video = 1; 

% Parameters that match our experiments with PAA
% Young's modulus
% Parameter.G = 230; 
% Parameter.G = 430; 
Parameter.G = 1500; 
% Poisson's ratio
Parameter.nu = 0.5; 
Parameter.E = 2*( 1+Parameter.nu )*Parameter.G; 

% Regularization parameter
lambda = Parameter.lambda; 

% Overlay the BFld image (1) or not (0)
flag_overlay = 1; 

% Frame rate of the output video
fr = Parameter.fr; 

% Filter size for the bilateral filter
sigma_s = 15; 
sigma_r = 10; 
Parameter.sigma_s = sigma_s; 
Parameter.sigma_r = sigma_r; 
% ======================================================================= %


% Input routes and file names
routeBFld = [route '\BFld\Images']; 
% routeLaser = [folder '\' fileName '\Beads_' num2str(chl) '\Images']; 
% routeIn = [folder '\' fileName '\Beads_' num2str(chl) '\DispField']; 
% PIVList = FUN_FileName(routeIn,'.mat','PIV_','No'); 

%% Read in displacement data
% Number of frames
N_TFM = size(Disp.U,3); 
% Grid size in PIV (pixel)
d = Disp.X(2)-Disp.X(1); 

% Time and length scales
lscl = Parameter.lscl(1)/1E6;   % Unit: m/pixel
tscl = Parameter.tscl;          % Unit: s/frame
Parameter.box_expt = lscl*d;    % Unit: m
Parameter.d = d;                % Unit: pixel


% Output routes for force data and images
routeOut = [route '\Force']; 
if ~exist( [routeOut '\Force_Images'],'file' )
    mkdir( routeOut ); 
    mkdir( [routeOut '\Force_Images'] ); 
end



%% Get information from the first frame
% load([routeIn '\' PIVList{1} '.mat']); 
X = Disp.X; 
Y = Disp.Y; 

% Make the dimensions of the matrices ODD numbers
if rem(length(X),2) == 0
    X(end) = []; 
end
if rem(length(Y),2) == 0
    Y(end) = []; 
end
D = [length(X) length(Y)]; 
N_box = [D(2) D(1)]; 


% [X_mesh, Y_mesh] = meshgrid(X,Y); 
% Mesh grid for the extrapolated force map
% [X_mesh,Y_mesh] = meshgrid( (1+1/d):1/d:D(1), (1+1/d)/d:1/d:D(2) ); 
% [X_mesh,Y_mesh] = meshgrid( 1:1:D(1), 1:1:D(2) ); 
% X_mesh = X_mesh*d*lscl; 
% Y_mesh = Y_mesh*d*lscl; 
xF = X*lscl*1E6; 
yF = Y*lscl*1E6; 

% Mesh grid for the vectors
% [x_mesh,y_mesh] = meshgrid(X,Y); 
% x_mesh = (x_mesh-1)*d+1; 
% y_mesh = (y_mesh-1)*d+1; 

% Parameter.lambda = lambda; 
Parameter.N_box = N_box; 
% Parameter.delta = delta; 



%% Calculate the corresponding Green's functions
[GF,~,~] = FUN_GenerateGreenFun( Parameter, ... 
    [folderG '\GreensFunctions'], lambda ); 

% k
% [K_x, K_y] = meshgrid(k_x, k_y); 


%% Get the force map for all frames
% Get information from the first BFld image
I0 = double( imread([routeBFld '\Img_0001.tif']) ); 
n = size(I0); 
I_mean =  mean( I0,'all' ); 
I_std = std( I0,0,'all' ); 
I_BFld = zeros( n(1),n(2),N_TFM ); 

% Mask: 0 - good pixels, 1 - bad pixels. 
tMask = Disp.ind_ex( 1:D(2),1:D(1),1:N_TFM ); 

Trac.Fx = zeros( D(2),D(1),N_TFM ); 
Trac.Fy = zeros( D(2),D(1),N_TFM ); 

% Go through all the frames
for i = 1:1:N_TFM
    % Read in the BFld image
    imgBFld = sprintf('Img_%04d',i);
    I_BFld(:,:,i) = imread([routeBFld '\' imgBFld '.tif']); 
    
    % Load in the information for force map
    U_m = Disp.U( 1:D(2),1:D(1),i ); 
    V_m = Disp.V( 1:D(2),1:D(1),i ); 
    [f_x,f_y] = disp2traction( U_m, V_m, GF ); 
    
    % Data
    Trac.X = X; 
    Trac.Y = Y; 
    Trac.Fx(:,:,i) = f_x; 
    Trac.Fy(:,:,i) = f_y;  
end

% Process 
I_BFld = medfilt3( double(I_BFld),[3 3 1] ); 
xI = [1,n(2)]*lscl*1E6; 
yI = [1,n(1)]*lscl*1E6;    
    
% Rescale the image for the transparency
% I_min = I_mean; 
% I_max = I_mean+3*I_std; 
I_min = I_mean-I_std; 
I_max = I_mean+2*I_std; 
I_resc = (I_BFld-I_min) ./ (I_max-I_min); 
I_resc(I_resc<0) = 0; 
I_resc(I_resc>1) = 1; 
    
    
    

%% Separate high and low frequency modes
%{
% Use gaussian filter to separate frequencies
win_s = 0.5; 
win_t = 50; 
win_end = 20; 

% Pad ends and gaussian filter to get the low frequency mode
% Displacement in x
Fx_lf = traction_lf( Trac.Fx,win_end,win_s,win_t ); 
% Displacement in y
Fy_lf = traction_lf( Trac.Fy,win_end,win_s,win_t ); 

% High frequency mode
Trac.Fx_hf = Trac.Fx-Fx_lf; 
Trac.Fy_hf = Trac.Fy-Fy_lf; 

% Remove the poor quality data
Trac.Fx(tMask == 1) = NaN; 
Trac.Fy(tMask == 1) = NaN; 
Trac.Fx_hf(tMask == 1) = NaN; 
Trac.Fy_hf(tMask == 1) = NaN; 
%}

disp( 'Filtering with bilateral filter' )

dI = size( Trac.Fx ); 
Fx_lf = zeros( dI ); 
Fy_lf = zeros( dI ); 

for i = 1:1:dI(1)
    for j = 1:1:dI(2)
        out = bilateral_filter_1d( Trac.Fx(i,j,:),sigma_s,sigma_r ); 
        out = bilateral_filter_1d( out,sigma_s,sigma_r ); 
        Fx_lf(i,j,:) = bilateral_filter_1d( out,sigma_s,sigma_r ); 

        out = bilateral_filter_1d( Trac.Fy(i,j,:),sigma_s,sigma_r ); 
        out = bilateral_filter_1d( out,sigma_s,sigma_r ); 
        Fy_lf(i,j,:) = bilateral_filter_1d( out,sigma_s,sigma_r ); 
    end
end

Trac.Fx_hf = Trac.Fx-Fx_lf; 
Trac.Fy_hf = Trac.Fy-Fy_lf; 

% Remove the poor quality data
Trac.Fx(tMask == 1) = NaN; 
Trac.Fy(tMask == 1) = NaN; 
Trac.Fx_hf(tMask == 1) = NaN; 
Trac.Fy_hf(tMask == 1) = NaN; 

% Save results
save([routeOut '\Trac.mat'],'Trac'); 



%% Make video
% Magnitude of traction
% F_mag = sqrt( Trac.Fx.^2+Trac.Fy.^2 ); 

% Traction field with overlay
videoName = 'Traction_overlay'; 
% v = VideoWriter( [routeOut '\' videoName '.avi'],'Uncompressed AVI' ); 
v = VideoWriter( [routeOut '\' videoName '.mp4'],'MPEG-4' ); 
v.FrameRate = fr; 
disp('Making traction video...')
open(v)

hd = figure('units','inch','position',[2,2,6,7.2]); 
for i = 1:1:N_TFM   
    cmap = gray(255); 
    RGB = ind2rgb( I_BFld(:,:,i),cmap );
    
    FUN_plot2DField( Trac.Fx(:,:,i),Trac.Fy(:,:,i), ... 
        xF,yF, 1, 'reverse', 'summer', 1, qScale ); 
    % imagesc( xF,yF,F_mag(:,:,i) ); 
    % imagesc( xI,yI,I_BFld(:,:,i) )
    hold on
    p2 = imshow(RGB, ... 
        'XData',[min(xI) max(xI)],...
        'YData',[min(yI) max(yI)]); 
    p2.AlphaData = I_resc(:,:,i).^1*0.6; 
    % p2.AlphaData = I_resc(:,:,i).^1*0.9; 
    hold off
    caxis( [0 3*std(sqrt(Trac.Fx.^2+Trac.Fy.^2),0,'all','omitnan')] )
    
    ax = gca; 
    ax.YDir = 'reverse'; 
    axis equal
    axis([0 max(xI) 0 max(yI)])
    % set(gcf,'color','w'); 
    % colormap('gray'); 
    colorbar;  
    
    %{
    % Save image
    imgName = sprintf('Img_%04d',i);
    % High resolution
    print('-dpng','-r200',[routeOut '\Force_Images\' imgName '.png'])
    %}
    
    % Save video
    ImgVideo = getframe(hd); 
    writeVideo(v, ImgVideo); 
    clf
end
close(v)

%{
% High frequency mode
videoName = 'Traction_hf'; 
v = VideoWriter( [routeOut '\' videoName '.avi'],'Uncompressed AVI' ); 
v.FrameRate = fr; 
disp('Making traction video...')
open(v)

hd = figure('units','inch','position',[2,2,5,6]); 
for i = 1:1:N_TFM       
    FUN_plot2DField( Trac.Fx_hf(:,:,i),Trac.Fy_hf(:,:,i), ... 
        xF,yF, 1, 'reverse', 'parula', 1, qScale);
    caxis([0 50])
    
    ax = gca; 
    ax.YDir = 'reverse'; 
    axis equal
    axis([0 max(xI) 0 max(yI)])
    set(gcf,'color','w'); 
    
    % Save video
    ImgVideo = getframe(hd); 
    writeVideo(v, ImgVideo); 
    clf
end
close(v)


% Total displacement
videoName = 'Traction_lf'; 
v = VideoWriter( [routeOut '\' videoName '.avi'],'Uncompressed AVI' ); 
v.FrameRate = fr; 
disp('Making traction video...')
open(v)

hd = figure('units','inch','position',[2,2,5,6]); 
for i = 1:1:N_TFM       
    FUN_plot2DField( Trac.Fx(:,:,i)-Trac.Fx_hf(:,:,i), ... 
        Trac.Fy(:,:,i)-Trac.Fy_hf(:,:,i), ... 
        xF,yF, 1, 'reverse', 'parula', 1, qScale); 
    caxis([0 50])
    
    ax = gca; 
    ax.YDir = 'reverse'; 
    axis equal
    axis([0 max(xI) 0 max(yI)])
    set(gcf,'color','w'); 
    
    % Save video
    ImgVideo = getframe(hd); 
    writeVideo(v, ImgVideo); 
    clf
end
close(v)
%}


% Save parameters
save([route '\Parameters.mat'],'Parameter'); 


% ==========================
end
% ==========================



%% Local functions
function F_lf = traction_lf( F,win_end,win_s,win_t )
    F_first = mean( F(:,:,1:win_end),3 ); 
    F_end = mean( F(:,:,end-win_end+1:end),3 ); 
    F = cat( 3,F_first,F,F_end ); 
    F_lf = imgaussfilt3( F,[win_s win_s win_t], ... 
        'FilterSize',[1 1 2*ceil(2*win_t)+1] );
    F_lf(:,:,[1 end]) = []; 
end 




