function Vectors = PIV_replaceOutliersBySecondPeak(PIVParams, Vectors)


    % read data
    x = Vectors.x;
    y = Vectors.y;
    dx = Vectors.dx;
    dy = Vectors.dy;
    dx2 = Vectors.dx2;
    dy2 = Vectors.dy2;
    validationFlag = Vectors.validationFlag;
    
    % first reshape everything in matrices
    % get dimensions of vector field
    i = find(y>y(1),1,'first');
    width = i-1;
    height = length(x) / width;
    % reshape vectors into matrices:
    dx = reshape(dx,width,height)';
    dy = reshape(dy,width,height)';
    validationFlag = reshape(validationFlag,width,height)';
    
    % find positions that need to be replaced
    [m, n] = find(validationFlag==1);
    % make a copy with zero's on the invalid positions
    dxOld = dx;
    dxOld(m,n) = 0;
    dyOld = dy;
    dyOld(m,n) = 0;
    % make a matrix with only ones on the positions of good vectors
    counter = ones(height,width);
    counter(m,n) = 0;
    % put a border of zeros around all these matrices, which makes handling
    % vectors on the edge much easier
    temp = zeros(height+2,width+2);
    temp(2:end-1,2:end-1) = dxOld;
    dxOld = temp;
    temp = zeros(height+2,width+2);
    temp(2:end-1,2:end-1) = dyOld;
    dyOld = temp;
    temp = zeros(height+2,width+2);
    temp(2:end-1,2:end-1) = counter;
    counter = temp;
    if ~isempty(m)
        for i=1:length(m)
            M = m(i);
            N = n(i);
            % we have to add 1 because we put a border of zeros around it:
            M = M + 1;
            N = N +1;
            % calculate average of nearest neighbours (= simplest form of
            % interpolation)
            % anywhere else in the image: 4 neighbours
            dxSum = (dxOld(M-1,N) + dxOld(M+1,N) + dxOld(M,N-1) + dxOld(M,N+1)) / 4;
            dySum = (dyOld(M-1,N) + dyOld(M+1,N) + dyOld(M,N-1) + dyOld(M,N+1)) / 4;
            numValidVectors = (counter(M-1,N) + counter(M+1,N) + counter(M,N-1) + counter(M,N+1)) / 4;
            % subtract 1 because we are dealing with the original matrices
            % again:
            M = M -1;
            N = N -1;
            if numValidVectors==0
                % not enough valid vectors around to interpolate
                dx(M,N) = 0;
                dy(M,N) = 0;
                % change flag to indicate that the vector still is wrong:
                validationFlag(M,N) = 1;
                % indicate that we're not done yet:
            else
                dx(M,N) = dxSum / numValidVectors;
                dy(M,N) = dySum / numValidVectors;
                % change flag to indicate that the vector was replaced:
                validationFlag(M,N) = 2;
            end
        end
    end
    
    % reshape matrices back into vectors:
    N = width * height;
    Vectors.dx = reshape(dx',1,N);
    Vectors.dy = reshape(dy',1,N);
    Vectors.validationFlag = reshape(validationFlag',1,N);
end