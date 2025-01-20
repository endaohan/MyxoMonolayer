%*************************************************************************%
% Read in images from the confocal microscopes in the format of .nd2
% Need the package bfmatlab
% Written by Endao Han, 2020/4/9, Version 6.
% Memory saving version. Fix the problem when a time series is saved as
% multiple series of images. 
%*************************************************************************%
% clear

%{
% ========================================================================%
%% Set up basic parameters
% Folder with data
% From computer
% folder = 'D:\DATA_Confocal\20190215'; 
% From removable drive
folder = 'F:\DATA_Confocal\20190327'; 
% folder = 'F:\DATA_Confocal\20190317'; 
% folder = 'G:\MyxoProject\DATA_Confocal\20190327'; 

% File name
fileList = FUN_FileName(folder,'.nd2','PAA','No'); 
N_file = length(fileList); 
fileName = fileList{7}; 
% fileName = fileList{6}; 

% Choose mode of data organization
% The data will have multiple colors. Each is one element in a cell. 
% Two modes: Z-Stack (0) and PIV (or time series) (1). 
% If Z-Stack is chosen, for each color, a cell array is created. Each
% element is one time step, and in each time step there is a 3D array whose
% three axes are x, y, and z. 
% If PIV is chosen, for each color, a cell array is created. Each element
% is one z level, and at each z there is a 3D array whose three axes are x,
% y, and time. 
flag_Mode = 0; 

% Save .mat data? Yes - 1; No - 0. 
flag_saveMat = 0; 

% ========================================================================%
%}

function [dataDim, dImg, DATA] = TFM_ReadImages(folder, fileName, Parameter)

if nargin == 3
    flag_Mode = Parameter.flag_Mode; 
    flag_saveMat = Parameter.flag_saveMat; 
elseif nargin == 2
    flag_Mode = 'TimeSeries'; 
    flag_saveMat = 0; 
end

% Channel names
ChName = {'BFld','Beads 1','Beads 2','Beads 3'}; 

%% Read in all the images and separate them based on the channels

% Read the .nd2 file
display(fileName)
data = bfopen([folder '\' fileName '.nd2']); 

% This should be the regular case. 
if size(data,1) == 1
    % Get information on color channel, z-stack, and time steps. 
    [nZ, nCh, N_img] = FUN_readZCT( data{1,1}{1,2} ); 
    dataDim = [nZ, nCh, N_img]; 

    % Decide what format to use to organize the data
    if N_img == 1 
        flag_Mode = 'ZStack'; 
        disp('The data is a Z-Stack. ')
    elseif nZ == 1
        flag_Mode = 'TimeSeries'; 
        disp('The data is a Time-series. ')
    else
        disp('The data is a Time-series of Z-Stack images. ')
    end

    % Get information from the first image
    dImg = size(data{1,1}{1,1}); 


    % Time series of Z-Stack
    if strcmp( flag_Mode, 'ZStack' )
        % Define an image sequence matrix for each channel
        DATA = cell(nCh,2);  

        % Split data, look at one channel and one slice at a time
        for i = 1:1:nCh
            display(['Now working on Channel ' num2str(i)])
            DATA{i,1} = cell(N_img,2);
            for j = 1:1:N_img
                DATA{i,1}{j,1} = zeros(dImg(1),dImg(2),nZ); 
                for k = 1:1:nZ
                    ind = (j-1)*nCh*nZ+(k-1)*nCh+i; 
                    DATA{i,1}{j,1}(:,:,k) = data{1,1}{ind,1}; 
                end
                DATA{i,1}{j,2} = data{1,1}{(j-1)*nCh*nZ+i,2}; 
            end
            DATA{i,2} = ChName{i}; 
        end


    % Z-Stacks of time series
    elseif strcmp( flag_Mode, 'TimeSeries' )

        % Define an image sequence matrix for each channel
        DATA = cell(nCh,2); 

        for i = 1:1:nCh
            display(['Now working on Channel ' num2str(i)])
            DATA{i,1} = cell(nZ,2); 
            for j = 1:1:nZ
                DATA{i,1}{j,1} = zeros(dImg(1),dImg(2),N_img); 
                for k = 1:1:N_img
                    ind = (k-1)*nCh*nZ+(j-1)*nCh+i; 
                    DATA{i,1}{j,1}(:,:,k) = data{1,1}{ind,1}; 
                end 
                DATA{i,1}{j,2} = data{1,1}{(j-1)*nCh+i,2}; 
            end
            DATA{i,2} = ChName{i}; 
        end

    else
        disp('EH ERROR: flag_Ch is wrong!!')
        return
    end
    
    
% This is when the .nd2 file is not saved correctly. 
% For now it only appeared for time-series data. 
else
    % Get information on color channel, z-stack, and time steps. 
    nZ = 1; 
    nCh = size(data{1,1},1); 
    N_img = size(data,1); 
    dataDim = [nZ, nCh, N_img]; 

    % Decide what format to use to organize the data
    disp('The data is a Time-series. ')

    % Get information from the first image
    dImg = size( data{1,1}{1,1} ); 
    % Define an image sequence matrix for each channel
    DATA = cell( nCh, 2 ); 

    for i = 1:1:nCh
        display(['Now working on Channel ' num2str(i)])
        DATA{i,1} = cell(nZ,2); 
        DATA{i,1}{1,1} = zeros(dImg(1),dImg(2),N_img); 
        for k = 1:1:N_img
            DATA{i,1}{1,1}(:,:,k) = data{k,1}{i,1}; 
        end 
        DATA{i,1}{1,2} = data{k,1}{i,2}; 
        DATA{i,2} = ChName{i}; 
    end
end


%% Save data as .mat file
if flag_saveMat == 1
    save([folder '\' fileName '.mat'],'flag_Mode','dataDim','DATA'); 
end
    


end




