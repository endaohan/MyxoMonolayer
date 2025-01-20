%*************************************************************************%
% TFM_main
% The main function for analyzing TFM results. 
% Written by Endao Han, V18, 2023/9/1
%*************************************************************************%
function TFM_main( folderIn, fileName )
% clear
% close all
% Data file name
% fileName = 'ExampleData'; 
% Route to the data
% folderIn = pwd; 

%% PART 1: Read in data
path_data = [folderIn '\' fileName]; 
% Load parameters
load( [path_data '\Parameters.mat'],'Parameter' )
% Path to subfunctions
addpath( [folderIn '\Function'],[folderIn '\track_ED'],[folderIn '\TFM_PIV'] )

% Run PIV analysis (1) or not (0)
flag_PIV = 1; 



%% PART 2: Run PIV        
if flag_PIV == 1
    disp( ['Running PIV on ' fileName] )

    % Create folders for PIV
    chl = Parameter.chl_PIV; 
    for j = chl
        % For TFM PIV (get displacement field)
        if ~exist([folderIn '\' fileName '\Beads_' num2str(j) '\TFM_Data'],'file')
            mkdir([folderIn '\' fileName '\Beads_' num2str(j) '\TFM_Images'])
            mkdir([folderIn '\' fileName '\Beads_' num2str(j) '\TFM_Data'])
            for k = 1:1:Parameter.nRef
                mkdir([folderIn '\' fileName '\Beads_' num2str(j) ... 
                    '\TFM_Data\TFM_Data_Ref' num2str(Parameter.iRef(k)) ])
            end
            mkdir([folderIn '\' fileName '\Beads_' num2str(j) '\TFM_Images_Col'])
        end
    end


    % Save the data obtained using every reference frame, but only save
    % PIV images when some reference frames are used to save time.  
    % When these are ref. frames, show the PIV images
    ind_Show = ((0:2:10)+1)*Parameter.tOverlap+1; 

    %% =============================================================
    % Write Parameter.txt files for PIV
    % Write the text files. 
    for j = chl
        routePIV = [folderIn '\' fileName '\Beads_' num2str(j)]; 

        for k = 1:1:Parameter.nRef
            if ismember( Parameter.iRef(k),ind_Show )
                FUN_PIV_readWriteParameters_TFM( routePIV,Parameter,Parameter.iRef(k),1 ); 
            else
                FUN_PIV_readWriteParameters_TFM( routePIV,Parameter,Parameter.iRef(k),0 ); 
            end
        end
    end

    % Perform PIV
    for j = chl
        Parameter.refLR = round( mean(Parameter.bd_TW) ); 
        routePIV = [folderIn '\' fileName '\Beads_' num2str(j)]; 

        % For loop to run PIV
        iRef = Parameter.iRef; 
        % For parallel computing
        p = parpool(4); 
        parfor k = 1:1:Parameter.nRef
            % Run PIV
            TFM_PIV_ref( [routePIV '\parameters_TFM_Ref' ... 
                num2str(iRef(k)) '.txt'], iRef(k) )
        end
        delete(p)
    end
    %% =============================================================

    % PIV is done. 
    disp('PIV is done.')
    close all
end



%% Post-analysis and presentation
% Calculate mean displacement
disp('Average displacement field. ')
% Use median filter (1) or not (0)
Parameter.flag_med = 0; 

chl = Parameter.chl_PIV; 
routeDrift = [folderIn '\' fileName '\Beads_' num2str(chl)]; 
% Get the displacement vector field
TFM_disp_new( routeDrift, Parameter, Parameter.flag_med ); 

% Find the shear modulus of the substrate
Parameter.G = 430; % Shear modulus is 430 Pa 
Parameter.lambda = 1/1E9; % Regularization parameter
TFM_disp2force_new( folderIn, fileName, Parameter, chl ); 

% Calculate the director field and find defects
BFld_finddefect( folderIn,fileName,0 ); 

% Make videos
% In the demo, only keep the high frequency component of traction. 
video_combined( folderIn,fileName,1 );  

% TFM is done. 
disp('TFM is done.')
close all
    

% ===== %
end
% ===== %





