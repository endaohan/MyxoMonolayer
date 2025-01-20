%*************************************************************************%
% File information
% Get file information in a folder
% Written by Endao Han, V5, 2020/4/9
% Now FUN_fileInfo only reads in information, but does not write. 
%*************************************************************************%

function [N_file, fileInfo] = FUN_fileInfo(excelName)

%{
clear
% For debugging
folder = 'D:\DATA_Confocal\20190515'; 
date = '20190515'; 
fileNameExample = 'PAA'; 
%}

% Save a .mat file? 1 - save; 0 - do not save.  
flagSave = 0; 

%% Get file information from the .xlsx file
% Read in the .xlsx file for information of the experiments
% recInfo: information recorded in the .xlsx file
recInfo = readcell(excelName); 

% Get rid of the "missing" boxes
for i = 1:1:size(recInfo,1)
    for j = 1:1:size(recInfo,2)
        if ismissing(recInfo{i,j})
            recInfo{i,j} = [];
        end
    end
end

% Number of files. The first row is column title, so exclude that. 
N_file = size(recInfo,1)-1; 

% Create fileInfo cell to store information about the files
fileInfo = recInfo; 


%% Save "fileInfo"
if flagSave == 1
    if ~exist([folder '\fileInfo.mat'],'file')
        % save([folder '\fileInfo.mat'],fileInfo); 
    end
end


end

