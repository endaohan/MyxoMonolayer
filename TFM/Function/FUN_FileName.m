% ======================================================================= %
% FUN_FileName: route is the folder where the files are; fileType is the
% extension of the files, e.g. txt, csv, etc; fileNameExample gives a
% format of the names that we want to look at.
% withExtension decides if the file name output include the extension name:
% Yes: with extension; No: without extension. 
% All inputs should be strings of characters. 
% Version 2. 2018/09/12
% ======================================================================= %

function file_name = FUN_FileName(route,fileType,fileNameExample,withExtension)

if strcmp(fileType,'folder')
    file_list = dir([route '/*']); 
else
% Find all the files with a certain file type in that folder
    file_list = dir([route '/*' fileType]);
end
name_list = cell(length(file_list),1);
ind_list = zeros(length(file_list),1);
indFolder = zeros(length(file_list),1);
indExample = zeros(length(file_list),1);

% Make sure that the files are 'PhiXX_X', which are how the files are named
for j = 1:1:length(file_list)
    name_list{j} = file_list(j).name;
    if strcmp(fileType,'folder')
        indFolder(j) = file_list(j).isdir;
        indExample(j) = ~isempty(strfind(name_list{j,1},fileNameExample));
        ind_list(j) = indFolder(j) && indExample(j);
    else
        ind_list(j) = ~isempty(strfind(name_list{j,1},fileNameExample));
    end
end
file_name = name_list(ind_list == 1); 

if ~strcmp(fileType,'folder')
    if strcmp(withExtension, 'No')
        for i = 1:1:length(file_name)
            name_temp = file_name{i};
            ind_dot = find(name_temp == char(46),1,'last'); 
            name_temp(ind_dot:end) = [];
            file_name{i} = name_temp;
        end
    end
end


