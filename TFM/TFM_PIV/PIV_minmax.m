function imgreturn = PIV_minmax(img,kernel)

% img = image to be processed
% kernel =  size of the grid which is min-maxed

IntImage = false;
if isinteger(img)
    IntImage = true;
    img = double(img);
end

if nargin<2
    kernel=31;
end

I_min = ones(size(img,1),size(img,2));
I_max = ones(size(img,1),size(img,2));

for x = ceil(kernel/2):size(img,2)-floor(kernel/2)        % loop over x with step 1
    %for y = ceil(kernel/2):kernel:size(img,1)-ceil(kernel/2)    % and over y
    for y = ceil(kernel/2):size(img,1)-floor(kernel/2)    % and over y
        %m = max(max(img(y-floor(kernel/2):y+floor(kernel/2),x-floor(kernel/2):x+floor(kernel/2))));     % m = maximum of this grid
        I = img(y-floor(kernel/2):y+floor(kernel/2),x-floor(kernel/2):x+floor(kernel/2));
        %get darkest pixel
        I_min(y,x) = min(min(I));
        %get brightest pixel
        I_max(y,x) = max(max(I));
    end
end

% smoothing needed!!!!!!!!!
% I_min = smooth(I_min)
% I_max = smooth(I_max)

% Smoothing filter
% h = ones(floor(kernel),floor(kernel)) / (floor(kernel))^2;
h = fspecial('gaussian');
I_min = imfilter(I_min,h);
I_max = imfilter(I_max,h);

imgreturn = img - I_min;
imgreturn = 255 * imgreturn ./ (I_max - I_min);

% imgreturn = img ./ (I_max - I_min);
if IntImage
%     imgreturn = img - I_min;
%     imgreturn = 255 * imgreturn ./ (I_max - I_min);
    imgreturn = uint8(imgreturn);
end