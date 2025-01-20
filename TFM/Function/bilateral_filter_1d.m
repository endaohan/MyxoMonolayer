% Test the 1d bilateral filter based on Gregory Kalliatakis
% Written by Endao Han, V1, 2/2/2021

function output = bilateral_filter_1d( I0,sigma_s,sigma_r )
%% BILATERAL FILTERING ON GRAYSCALE IMAGES

% Illustrates the use of bilateral filtering on grayscale images.
% Bilateral filtering is a technique to smooth images while preserving edges, 
% by means of a nonlinear combination of nearby
% image values. The method is noniterative, local, and simple.
% The bilateral filter is controlled by two parameters: ?s (spatial parameter) and ?r (range parameter)

% parameters:
% inputImg    = input image
% sigma_s  = domain parameter for spatial kernel
% sigma_r  = range parmeter for intensity kernel
% noise     = Gaussian noise intensity

% Code Developed by :
% Gregory Kalliatakis (December 2014) - gkalliatakis@yahoo.gr
% Masters in Computer Vision
% University of Burgundy, France
% Advanced Image Analysis - Homework Project on Bilateral Filtering

%% Example of usage
% im=imread('images/lena.jpg');
% out=bilateralGrayscale(im,3,0.2,0.1);

%% Read in example data
%{
clear
close all

routeIn = 'D:\DATA_Confocal\20210126\PAA430Pa_Beads100nm_noCover_Coat_TimeSeries_Depth_60x_1x_1\Force'; 
load( [routeIn '\Trac.mat'] ); 

% I0 = Trac.Fx(200,50,:); 
I0 = Trac.Fx(100,50,:); 
% figure; plot(I0)
%}



%% Set parameters
% sigma_s = 30; 
% sigma_r = 30;
% image=double(image)/255;
% window_size=9; 
window_size = 2*ceil(2*sigma_s)+1; 

if length(size(I0)) == 3
    [row, col, d]=size(I0); 
    if (row~=1)||(col~=1)
        disp('EROR: bilateral_filter_1d input not a 1d array!')
        return
    end
    I0 = reshape( I0,[d 1] ); 
else
    I0 = reshape( I0,[length(I0) 1] ); 
end
    
noisy_image = [ ones(1,window_size)*I0(1) I0.' ones(1,window_size)*I0(end) ].'; 

%%
% transforms the domain specified by vectors from -window_size to window_size into arrays X and Y 
x = meshgrid( -window_size:window_size,1 ).';

%% Domain filter 
%The weights depend on the spatial distance (to the center pixel x) only; therefore, it is calculated once and saved.
domain_filter=exp(-x.^2/(2*sigma_s^2));  
% domain_filter = x.^0/length(x);  

%% Repeat for all pixels
r = size(noisy_image,1);
output = zeros(size(noisy_image));

for i=1:r
    %Adjusting the window size
    imin=max(i-window_size,1);
    imax=min(i+window_size,r);
    I = noisy_image(imin:imax,1);

    range_filter=exp(-double(I-noisy_image(i,1)).^2/(2*sigma_r^2)); % Range filter

    %Taking the product of the range and domain filter.The combination is refered to as Bilater Filter
    BilateralFilter=range_filter.*domain_filter( (imin:imax)-i+window_size+1,1 );

    Fnorm=sum(BilateralFilter(:));
    output(i,1)=sum(sum(BilateralFilter.*double(I)))/Fnorm; %normalize the output
end

output( [1:window_size end-window_size+1:end] ) = []; 

%{
figure; 
plot(I0)
hold on
plot(output(window_size+1:end-window_size),'r-')
hold off
%}

end

