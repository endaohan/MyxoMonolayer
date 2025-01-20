function Vectors = PIV_correlate(PIVParams, Vectors, Im1, Im2)

%
% PIV_correlate - Correlate two images according to the parameters
% (This function is called by the function PIV)
%


step = PIVParams.GridSize - PIVParams.Overlap;
grid = PIVParams.GridSize; %This is a bit shorter than everytime PIVParams.GridSize :)

% get the size of the images
[height, width] = size(Im1);

% set the positions of the corners of our correlation windows:
if PIVParams.UseROI
    % we have a helper function which creates the correct coordinates for us:
    [x_points, y_points, ~] = PIV_ROI(PIVParams);
    %     % set points only inside the ROI (Region Of Interest)
    %     y_points = PIVParams.ROITopEdge:step:PIVParams.ROIBottomEdge-grid;
    %     x_points = PIVParams.ROILeftEdge:step:PIVParams.ROIRightEdge-grid;
else
    % use the complete image
    y_points = 1:step:height-grid+1;
    x_points = 1:step:width-grid+1;
end


%% The loop going through all the correlation windows
index = 1;
for y = y_points
    for x = x_points
        % pick one cell of the grid:
        y1 = y;
        y2 = y+grid-1;
        x1 = x;
        x2 = x+grid-1;
        ImA = Im1(y1:y2,x1:x2);
        % from the second image we take the same cell, with a
        % displacement if this is not the first pass (for the first
        % pass the displacements are defined as zero)
        y1 = y1 + round(Vectors.dy(index));
        y2 = y2 + round(Vectors.dy(index));
        x1 = x1 + round(Vectors.dx(index));
        x2 = x2 + round(Vectors.dx(index));
        %Check if patches are not outside the image. If so, just take
        %the original patch without shift.
        if y1<1 || y2>height
            y1 = y1 - round(Vectors.dy(index));
            y2 = y2 - round(Vectors.dy(index));
            Vectors.dy(index)=0;
        end
        if x1<1 || x2>width
            x1 = x1 - round(Vectors.dx(index));
            x2 = x2 - round(Vectors.dx(index));
            Vectors.dx(index)=0;
        end
        % Check for other weird things
        if isnan(x1) || isnan(x2)
            x1 = x;
            x2 = x + grid-1;
        end
        if isnan(y1) || isnan(y2)
            y1 = y;
            y2 = y + grid-1;
        end
        ImB = Im2(y1:y2,x1:x2);
        
        
        
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Correlation with XCORR method %%%%%%
        if strncmp(PIVParams.Method, 'XCORR', 5)
            %normalize images
            % The following 4 steps are needed to get a meaningful
            % value for the correlation coefficient (see "PIV a
            % practical guide", Second Edition, Springer)
            % step 1: compute mean and standard deviation of samples
            ImAMean = mean(ImA(:));
            ImBMean = mean(ImB(:));
            ImASTD = std(ImA(:));
            ImBSTD = std(ImB(:));
            % step 2: subtract mean from samples
            ImA = ImA - ImAMean;
            ImB = ImB - ImBMean;
            % step 3: calculate the cross-correlation of the samples
            if grid>24
                % for larger correlation windows, cross-correlation using
                % fft is faster
                CorrIm = fftxcorr2(ImB, ImA);
            else
                % for smaller correlation windows xcorr2 is faster
                CorrIm = xcorr2(ImB, ImA);
            end
            % step 4: divide by standard deviations
            CorrIm = CorrIm / (ImASTD*ImBSTD);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Correlation with NORMXCORR method %%%%%%
        if strncmp(PIVParams.Method, 'NORMXCORR', 9)
            % Calculate the normalized cross-correlation of the samples
            CorrIm = normxcorr2(ImB, ImA);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        
        %----------------------------------------------
        % Peak detection
        if PIVParams.InGridRefinement || ~PIVParams.UseSecondPeak
            % If we are still in the process of grid-refinement there is no
            % need to do multiple peak detection
            % single peak detection
            corrValue = max(CorrIm(:));
            [Y, X] = find(CorrIm==corrValue);
            if numel(X)
                X = X(1);
                Y = Y(1);
            else
                X = grid;
                Y = grid;
            end
            X2 = grid;
            Y2 = grid;
        else
            % multiple peak detection:
            [XX, YY, Values] = PIV_peakDetection(CorrIm);
            Y = YY(1);
            X = XX(1);
            corrValue = Values(1);
            if length(XX)>1
                X2 = XX(2);
                Y2 = YY(2);
            else
                X2 = grid;
                Y2 = grid;
            end
        end
        
        
        
        
        
        %------------------------------------------------
        % subpixel peak detection
        
        % subpixel displacement first peak
        if (Y<=1 || Y==(2*grid-1)) || (X<=1 || X==(2*grid-1))
            % The peak is on an edge pixel of the correlation
            % window, so we cannot calculate a subpixel
            % displacement.
            deltaY = 0;
            deltaX = 0;
        else
            [deltaY, deltaX] = PIV_subpixelpeak(CorrIm(Y-1:Y+1,X-1:X+1));
        end
        % displacement vector:
        dy = Y + deltaY - grid;
        dx = X + deltaX - grid;
        % if we don't find a number we might as well use 0
        if isnan(dy)
            dy=0;
        end
        if isnan(dx)
            dx=0;
        end
        
        % subpixel displacement second peak
        if Y2==0 || Y2==1 || Y2==(2*grid-1) || X2==0 || X2==1 || X2==(2*grid-1)
            % The peak is on an edge pixel of the correlation
            % window, so we cannot calculate a subpixel
            % displacement.
            deltaY = 0;
            deltaX = 0;
        else
            [deltaY, deltaX] = PIV_subpixelpeak(CorrIm(Y2-1:Y2+1,X2-1:X2+1));
        end
        % displacement vector:
        dy2 = Y2 + deltaY - grid;
        dx2 = X2 + deltaX - grid;
        % if we don't find a number we might as well use 0
        if isnan(dy2)
            dy2=0;
        end
        if isnan(dx2)
            dx2=0;
        end
        
        
        
        
        
        
        
        %--------------------------------------------------
        % save data
        % register current position:
        Vectors.y(index) = y + round(grid/2);
        Vectors.x(index) = x + round(grid/2);
        % register displacement vector:
        Vectors.dx(index) = round(Vectors.dx(index)) + dx;
        Vectors.dy(index) = round(Vectors.dy(index)) + dy;
        % register second peak displacement vector:
        Vectors.dx2(index) = round(Vectors.dx(index)) + dx2;
        Vectors.dy2(index) = round(Vectors.dy(index)) + dy2;
        % register correlation value:
        Vectors.correlationCoefficient(index) = corrValue;
        % increase index
        index = index + 1;
    end
end


    function C = fftxcorr2(A,B)
        % get sizes of the two matrices
        [mA,nA] = size(A);
        [mB,nB] = size(B);
        % the size of our output matrix
        m = mA + mB - 1;
        n = nA + nB - 1;
        % zero-pad the original matrices to match with the size of our
        % output matrix
        A1 = zeros(m,n);
        B1 = zeros(m,n);
        A1(1:mA,1:nA) = A;
        B1(mA:end,nA:end) = B;
        % do the cross correlation using fft transform
        C = ifft2(conj(fft2(A1)).*fft2(B1));
        % rotate the result
        C = rot90(C,2);
    end


end