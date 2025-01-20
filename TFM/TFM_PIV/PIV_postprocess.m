function Vectors = PIV_postprocess(PIVParams, Vectors)

%
% PIV_postprocess - post-process the PIV-vectors according to the
% parameters
% (This function is called by the function PIV)
%

%% Number of validationpasses
numValidationPasses = PIVParams.numValidationPasses;
PIVParams.validationPass = 1;
while PIVParams.validationPass<numValidationPasses+1
    
    %%% Mark outliers
    Vectors = PIV_validate(PIVParams, Vectors);
    
    %%% try replacing outliers by second peak
    if PIVParams.UseSecondPeak
        Vectors.dx(Vectors.validationFlag==1) = Vectors.dx2(Vectors.validationFlag==1);
        Vectors.dy(Vectors.validationFlag==1) = Vectors.dy2(Vectors.validationFlag==1);
        % Mark the vector that we used a second-order peak:
        Vectors.validationFlag(Vectors.validationFlag==1) = 3;
    end
    
    %%% again, mark the outliers
    Vectors = PIV_validate(PIVParams, Vectors);
    
    %%% replace remaining outliers using interpolation (if possible)
    if PIVParams.UseInterpolation
        Vectors = PIV_interpolateVectors(Vectors);
    end
    
    %%% Check if we have replaced all spurious data
    if isempty(find(Vectors.validationFlag==1,1))
        % this will get us out of the while loop:
        PIVParams.validationPass = numValidationPasses + 1;
    else
        % on to the next validation pass:
        PIVParams.validationPass = PIVParams.validationPass + 1;
    end
end

%% Data smoothing
if PIVParams.UseDataSmoothing
    % read data
    x = Vectors.x;
    y = Vectors.y;
    dx = Vectors.dx;
    dy = Vectors.dy;
    
    % first reshape everything in matrices
    % get dimensions of vector field
    width = find(y==y(1),1,'last');
%     width = i-1;
    height = length(x) / width;
    % reshape vectors into matrices:
    dx = reshape(dx,width,height)';
    dy = reshape(dy,width,height)';
    
    % create the gaussian smoothing filter of the right size:
    filterSize = PIVParams.SmoothingFilterSize;
    filterType = PIVParams.SmoothingFilterType;
    filter = fspecial(filterType,filterSize);
    % filter the data:
    dx = imfilter(dx,filter,'replicate');
    dy = imfilter(dy,filter,'replicate');
    
    % reshape matrices back into vectors:
    N = width * height;
    Vectors.dx = reshape(dx',1,N);
    Vectors.dy = reshape(dy',1,N);
end
    