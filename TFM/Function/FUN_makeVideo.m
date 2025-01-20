% Function to make a video
% Input: Direction to a folder with images. 
% Output: .mp4 file
% Written by Endao Han, V4, 2021/1/3. 
% Change the colormap for every channel to be the same. 
% ***********************************************************************

function FUN_makeVideo(videoIn,videoOut,videoName,fr,channelColor)

%{ 
% For debugging
clear

route = 'G:\DATA_Confocal\20190717\PAA230Pa_CTTYE_DpilA_noCover_TimeSeries_60x_1p5x_2'; 
videoIn = [route '\Beads_1\Images_Col'];
videoOut = route; 
videoName = 'video_Laser1'; 
fr = 20; 
channelColor = 'Beads_1'; 
%}            

%% Input parameters
% Input and output routes for the images
imgInfo.routeIn = videoIn; 
imgInfo.routeOut = videoOut; 
imgInfo.imgType = 'tif'; 
% imgInfo.imgType = 'png'; 
% Output route for the video
videoRouteOut = imgInfo.routeOut; 

% Image step
ImgSpacing = 1; 

%% Get folder information
% Find all the .txt files in a folder
fileList = dir([imgInfo.routeIn '/*.' imgInfo.imgType]);
N_img = length(fileList); 
N_video = floor(N_img/ImgSpacing); 

% Set colormap
%{
switch channelColor
    case 'BFld'
        cMap = gray(256); 
    case 'Beads_1'
        cMap = hot(256); 
    case 'Beads_2'
        cMap = cool(256); 
    case 'Beads_3'
        cMap = summer(256); 
    otherwise
        cMap = parula(256); 
end
%}
switch channelColor
    case 'BFld'
        cMap = gray(256); 
    case 'Beads_1'
        cMap = hot(256); 
    case 'Beads_2'
        cMap = hot(256); 
    case 'Beads_3'
        cMap = hot(256); 
    otherwise
        cMap = parula(256); 
end




%% OUTPUT
% Make a video
v = VideoWriter([videoRouteOut '\' videoName '.mp4'],'MPEG-4'); 
v.FrameRate = fr; 
open(v)
disp('Making video...')
for i = 1:1:N_video
    j = (i-1)*ImgSpacing + 1; 
    imgInfo.fileName = fileList(j).name; 
    ImgVideo = imread([imgInfo.routeIn '\' imgInfo.fileName],'tif'); 
    ImgRGB = ind2rgb(ImgVideo,cMap); 
    writeVideo(v, ImgRGB); 
end
close(v)

close all


