function Vectors = PIV_validate(PIVParams, Vectors)


%% Validation
% We use three different validation flags:
% 0: Valid vector
% 1: Invalid vector
% 2: Replaced vector (see replacement scheme)





%% The simplest test: maximum displacement
% This is a good one to do always just to get rid of ridiculous vectors
if PIVParams.UseMaxDisplacement
    displacement = sqrt((Vectors.dx).^2 + (Vectors.dy).^2);
    Vectors.validationFlag(displacement>PIVParams.MaxDisplacement) = 1;
end








%% This one only makes sense to do in the first validation pass:
% It is doubtful how much sense this test makes. A low correlation
% coefficient does not mean that the vector is false. Also, a spurious
% vector might very well have a large correlation coefficient.
if PIVParams.UseMinCorrelationCoefficient && PIVParams.validationPass==1
    corrValue = Vectors.correlationCoefficient;
    Vectors.validationFlag(corrValue<PIVParams.MinCorrelationCoefficient) = 1;
end







%% A good test: Normalized Median Test (see PIV-book)
% This is a very robust method. A bit more involved, but usually it's worth
% it.
if PIVParams.UseNormalizedMedianTest
    % read data
    x = Vectors.x;
    y = Vectors.y;
    dx = Vectors.dx;
    dy = Vectors.dy;
    validationFlag = Vectors.validationFlag;
    
    % first reshape everything in matrices
    % get dimensions of vector field
    width = find(y==y(1),1,'last');
%     width = i;
    height = length(x) / width;
    % reshape vectors into matrices:
    dx = reshape(dx,width,height)';
    dy = reshape(dy,width,height)';
    validationFlag = reshape(validationFlag,width,height)';
    threshold = PIVParams.MedianThreshold;
    
    % put a border of zeros around all these matrices, which makes handling
    % vectors on the edge much easier
    temp = zeros(height+2,width+2);
    temp(2:end-1,2:end-1) = dx;
    dx = temp;
    temp = zeros(height+2,width+2);
    temp(2:end-1,2:end-1) = dy;
    dy = temp;
    % a matrix with ones on the positions where we have vectors
    counter = ones(height,width);
    temp = zeros(height+2,width+2);
    temp(2:end-1,2:end-1) = counter;
    counter = temp;
    
    for m=1:height
        for n=1:width
            % we have to add 1 because we put a border of zeros around it:
            M = m + 1;
            N = n +1;
            % make lists of local dx's and dy's
            dx_list = dx(M-1:M+1,N-1:N+1);
            dx_list = dx_list(:);
            dy_list = dy(M-1:M+1,N-1:N+1);
            dy_list = dy_list(:);
            counter_list = counter(M-1:M+1,N-1:N+1);
            counter_list(2,2) = 0;
            counter_list = counter_list(:);
            % select the vectors that we want to use for our
            % comparison:
            dx_list = dx_list(counter_list==1);
            dy_list = dy_list(counter_list==1);
            % calculate median of dx and dy
            dx_med = median(dx_list);
            dy_med = median(dy_list);
            % calculate residuals
            dx_r = abs(dx_list - dx_med);
            dy_r = abs(dy_list - dy_med);
            % calculate median of residuals
            dx_r_med = median(dx_r);
            dy_r_med = median(dy_r);
            % epsilon_0
            epsilon_0 = 0.1;
            % do the test:
            x_test = abs(dx_med - dx(M,N)) / (dx_r_med + epsilon_0);
            y_test = abs(dy_med - dy(M,N)) / (dy_r_med + epsilon_0);
            if x_test > threshold || y_test > threshold
                % set vector to invalid (=1)
                validationFlag(m,n) = 1;
            end
        end
    end
    
    % reshape validation matrix back into vector:
    N = width * height;
    Vectors.validationFlag = reshape(validationFlag',1,N);
end