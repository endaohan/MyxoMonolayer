% main_PIV_Agar
% Process images and perform PIV to get cell velocity
% Based on main_PIV_PAA.m and main_velocity_compare
% 
% Written by Endao Han, V1, 9/1/2023
% 
function fun_velocity( fileName )
% For debugging
%{
clear
close all

fileName = 'ExampleData'; 
%}

%% Basic parameters
% Input and output routes
routeIn = [pwd '\' fileName]; 
routeOut = routeIn; 

addpath( 'Function' )

% Load parameters
lscl = 0.11;       % unit: micron
tscl = 15;       % unit: minute


%% Read in data
% Get basic parameters
fileList = dir( [routeIn '\BFld\*.tif'] ); 
N_img = size( fileList,1 ); 
I_0 = imread( [routeIn '\BFld\' fileList(1).name] ); 
dImg = [ size( I_0,1 ) size( I_0,2 ) ]; 
x = 1:1:dImg(2); 
y = 1:1:dImg(1); 

% Read in images
I_proc = nan( dImg(1),dImg(2),N_img ); 
I = nan( dImg(1),dImg(2),N_img ); 
for i = 1:1:N_img
    I_proc(:,:,i) = imread( [routeIn '\BFld\' fileList(i).name] ); 
    I(:,:,i) = imread( [routeIn '\Laser\' fileList(i).name] ); 
end
I = rescale(I); 


%% Optical flow of BFld images
X = 21:20:1001; 
Y = 21:20:1421; 

display( 'Calculating cell velocity using optical flow. ' )
[Vx_op,Vy_op] = cal_velocity_optFlow( I_proc,X,Y ); 
% Convert unit to micron/min
optFlow.X = X; 
optFlow.Y = Y; 
optFlow.Vx_op = Vx_op*lscl/(tscl/60); 
optFlow.Vy_op = Vy_op*lscl/(tscl/60); 
save( [routeOut '\Data_optFlow.mat'],'optFlow' ); 


% ==================
end
% ==================


