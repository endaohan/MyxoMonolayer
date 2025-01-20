function [PIVParams, Vectors] = PIV_gridRefinement(PIVParams, Im1, Im2)

% list of refinement steps:
gridSizes = [1 2 4 6 8 12 16 24 32 48 64 96 128 192 256 384];

% tell everyone that we are inside the grid refinement scheme:
PIVParams.InGridRefinement = true;

%% override some user settings that we never/always use in grid refinement
% save user settings to recover after grif refinement:
Overlap = PIVParams.Overlap;
UseInterpolation = PIVParams.UseInterpolation;
UseDataSmoothing = PIVParams.UseDataSmoothing;
SmoothingFilterType = PIVParams.SmoothingFilterType;
% and set the values we use for grid refinement:
PIVParams.Overlap = 0;
PIVParams.UseDataSmoothing = true;
PIVParams.SmoothingFilterType = 'average';
PIVParams.UseInterpolation = true;


%% save target grid size:
GridSize = PIVParams.GridSize;
PIVParams.GridSize = PIVParams.InitialGridSize;

%% we have to redo the pre-allocation of matrices
Vectors = PIV_generateStructure(PIVParams, Im1);

%% Do the refinement steps
while PIVParams.InGridRefinement
    display(['Grid: ' num2str(PIVParams.GridSize)])
    % correlate images
    Vectors = PIV_correlate(PIVParams, Vectors, Im1, Im2);
    % postprocess Vectors
    Vectors = PIV_postprocess(PIVParams, Vectors);
    % grid size for our next refinement step
    idx = find(gridSizes<PIVParams.GridSize,1,'last');
    NewGridSize = gridSizes(idx);
%     NewGridSize = uint16(2*round(PIVParams.GridSize / 4));
    % check if we need to further refine
    if NewGridSize <= GridSize
        % no: tell everyone that we are not inside the grid refinement
        % scheme anymore and get out of the while loop
        PIVParams.InGridRefinement = false;
    else
        % refine grid
        [PIVParams, Vectors] = PIV_refineGrid(PIVParams, Vectors, Im1, NewGridSize);
    end
end

%% set back the replacement and filter settings:
PIVParams.Overlap = Overlap;
PIVParams.UseInterpolation = UseInterpolation;
PIVParams.UseDataSmoothing = UseDataSmoothing;
PIVParams.SmoothingFilterType = SmoothingFilterType;
% set grid size to the target grid size:
[PIVParams, Vectors] = PIV_refineGrid(PIVParams, Vectors, Im1, GridSize);