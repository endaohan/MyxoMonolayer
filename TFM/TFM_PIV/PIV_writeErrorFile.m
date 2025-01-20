function PIV_writeErrorFile(PIVParams)

% format: each line consists of a string (%s) and newline command (\r\n)
format = '%s\r\n';

% Data filename for current frame:
FileName = [PIVParams.DataDir '_ERROR_frame' num2str(PIVParams.imageNumber,'%06u') '.txt'];

% open file
fid = fopen(FileName,'a');

% write info to file
fprintf(fid, format, 'PIV error log');
fprintf(fid, format, ['generated on: ' datestr(now,'mmmm dd, yyyy HH:MM:SS AM')]);
fprintf(fid, format, ['Frame number: ' num2str(PIVParams.imageNumber)]);
fprintf(fid, format, ['Message: ' PIVParams.AbortMessage]);
fprintf(fid, format, 'end of error message');
fprintf(fid, format, '********************');
fprintf(fid, '\r\n');
% close file
fclose(fid);