function [PIVParams, Im1, Im2] = PIV_prepocess_TFM(PIVParams, Im1, Im2)

%
% PIV_prepocess - preprocessing function for PIV
% (This function is called by the function PIV)
%

%% Mask image
if exist([PIVParams.Directory 'mask.png'],'file')
    mask = imread([PIVParams.Directory 'mask.png']);
    Im1 = Im1 .* uint8(logical(mask));
    Im2 = Im2 .* uint8(logical(mask));
%{
% This part was used temporarily for analyzing Jiepan's data
else
    mask = Im2-10;
    mask(mask<0) = 0; 
    Im1 = Im1 .* uint8(logical(mask));
    Im2 = Im2 .* uint8(logical(mask));
%}
end

%% Crop image

[PIVParams, Im1] = PIV_crop(PIVParams, Im1);
[PIVParams, Im2] = PIV_crop(PIVParams, Im2);

%% Convert the images to double precision
Im1 = double(Im1);
Im2 = double(Im2);

%% Image filters

%%% MinMax Filter
if PIVParams.UseMinMax
    Im1 = PIV_minmax(Im1,PIVParams.MinMaxKernel);
    Im2 = PIV_minmax(Im2,PIVParams.MinMaxKernel);
end

%%% Unsharp mask (sharpen image)
if PIVParams.UseUnsharpMask
    f = fspecial('gaussian');
    % create unsharp mask:
    Mask1 = imfilter(Im1,f,'replicate');
    Mask2 = imfilter(Im2,f,'replicate');
    % subtract mask from images:
    Im1 = Im1 - PIVParams.UnsharpMaskAmount*Mask1;
    Im2 = Im2 - PIVParams.UnsharpMaskAmount*Mask2;
%     Im1 = imsharpen(Im1,'Amount',PIVParams.UnsharpMaskAmount);
%     Im2 = imsharpen(Im2,'Amount',PIVParams.UnsharpMaskAmount);
end

%%% Gaussian blur
if PIVParams.UseGaussianBlur
    f = fspecial('gaussian',5,PIVParams.GaussianBlurRadius);
    Im1 = imfilter(Im1,f,'replicate');
    Im2 = imfilter(Im2,f,'replicate');
end

% figure
% I = Im1;
% I = I - min(I(:));
% I = I / max(I(:));
% I = uint8(I * 255);
% imshow(I)
% pause
% close


