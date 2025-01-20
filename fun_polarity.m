% main_polarity_Agar
% Main function of polarity measurement of cell layer
% Based on main_polarity_Agar.m
% 
% Written by Endao Han, V5, 8/15/2022
% 
function fun_polarity( fileName )
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
mkdir( [routeOut '\Label'] )
mkdir( [routeOut '\Polarity'] )

addpath( 'Function' )

% Load parameters
lscl = 0.11;       % unit: micron
tscl = 15;       % unit: minute

% Polarity of each cell: mglA - 1, mglB - -1
p_cell = -1; 

% Threshold for binarize image
I_thresh = 0.8; 



%% Read in data
% Get basic parameters
fileList = dir( [routeOut '\BFld\*.tif'] ); 
N_img = size( fileList,1 )-1; 
I_0 = imread( [routeOut '\BFld\' fileList(1).name] ); 
dImg = [ size( I_0,1 ) size( I_0,2 ) ]; 
x = 1:1:dImg(2); 
y = 1:1:dImg(1); 

% Read in images
I_proc = nan( dImg(1),dImg(2),N_img ); 
I = nan( dImg(1),dImg(2),N_img ); 
for i = 1:1:N_img
    I_proc(:,:,i) = imread( [routeOut '\BFld\' fileList(i).name] ); 
    I(:,:,i) = imread( [routeOut '\Laser\' fileList(i).name] ); 
end
I = rescale(I); 

% Read in optical flow data
load( [routeOut '\Data_optFlow.mat'],'optFlow' )
X = optFlow.X; 
Y = optFlow.Y; 
Vx_op = optFlow.Vx_op; 
Vy_op = optFlow.Vy_op; 



%% Measure director field and find defects
[dirf_all,ADefs,pdList,ndList,lays] = PV_finddefect( I_proc,routeOut ); 


%% Find centers of fluorescent labels
I_bw = imbinarize( I,I_thresh ); 
I_bw = bwareaopen( I_bw,5 ); 
ct = cell( N_img,2 ); 
for j = 1:1:N_img
    stats = regionprops( I_bw(:,:,j),'Area','Centroid' ); 
    a = struct2cell( stats );
    ct{j,1} = cell2mat( a(2,:).' ); 
    ct{j,2} = cell2mat( a(1,:).' );
end

figure; 
subplot(1,2,1)
imagesc(I(:,:,10))
hold on
plot( ct{10,1}(:,1),ct{10,1}(:,2),'r.' )
hold off
subplot(1,2,2)
imagesc(I_bw(:,:,10))



%% Plot overlaid images
hd_overlay = figure( 'Unit','Inch','Position',[2 2 4.5 6] ); 
for i = 1:1:N_img
    imagesc( I_proc(:,:,i) )
    colormap('gray')
    hold on
    % Fluorescent label
    plot( ct{i,1}(:,1),ct{i,1}(:,2),'g.' )
    % Defects
    for j = 1:1:size( pdList,1 )
        if ismember( i,pdList{j,1}.i )
            ip = find( pdList{j,1}.i == i ); 
            plot_defect_symbol( [pdList{j,1}.x(ip) pdList{j,1}.y(ip)], ... 
                pdList{j,1}.d(ip,:)*50,1 )
            text( pdList{j,1}.x(ip)+20,pdList{j,1}.y(ip)+20, ... 
                num2str( pdList{j,1}.id ),'Color','r' )
        end
    end
    for j = 1:1:size( ndList,1 )
        if ismember( i,ndList{j,1}.i )
            ip = find( ndList{j,1}.i == i ); 
            plot_defect_symbol( [ndList{j,1}.x(ip) ndList{j,1}.y(ip)], ... 
                ndList{j,1}.d(ip,:)*50,-1 )
            text( ndList{j,1}.x(ip)+20,ndList{j,1}.y(ip)+20, ... 
                num2str( ndList{j,1}.id ),'Color','b' )
        end
    end
    hold off
    axis equal
    axis( [0 dImg(2) 0 dImg(1)] )
    
    name_tiff = sprintf( 'Img_%04d',i ); 
    print( hd_overlay,[routeOut '\Label\' name_tiff '.png'],'-dpng','-r200'); 
    clf
end
FUN_ImgToVideo( [routeOut '\Label'], 'png', routeOut,'Label_video',13,'mp4' );



%% Calculate polarity
% Output: Plty
% Column: 1 - x, 2 - y, 3 - polarity x, 4 - polarity y, 5 - asymmetry, 6 - area
Plty = cell( N_img,1 ); 
w = 5; 
for i = 1:1:N_img
    if rem(i,10) == 0
        disp( [num2str(i) '/' num2str(N_img)] )
    end
    ct_1 = ct{i,1}; 
    ct_r = round( ct_1 ); 
    dirf = dirf_all(:,:,i);  
    p_asym = nan( size( ct_1,1 ),1 ); 
    p_sign = nan( size( ct_1,1 ),1 ); 
    p_vec = nan( size( ct_1,1 ),2 ); 
    I_pad = padarray( I_proc(:,:,i),[2*w 2*w] ); 
    % tic
    for j = 1:1:size( ct_1,1 )
        % [xq,yq,sfr] = rotate_scalar_field( x,y,I_proc(:,:,i), ... 
        %     w,ct_1(j,:),dirf( ct_r(j,2),ct_r(j,1) ) ); 
        [xq,yq,sfr] = rotate_scalar_field( 1:1:4*w+1,1:1:4*w+1, ... 
            I_pad( ct_r(j,2):ct_r(j,2)+w*4,ct_r(j,1):ct_r(j,1)+w*4 ), ... 
            w,[2*w+1 2*w+1],dirf( ct_r(j,2),ct_r(j,1) ) ); 
        b_sym_r = mean( sfr(w+1-2:w+1+2,w+2:end-1),'all','omitnan' ); 
        b_sym_l = mean( sfr(w+1-2:w+1+2,2:w),'all','omitnan' ); 
        p_asym(j,1) = b_sym_l-b_sym_r; 
        
        theta_temp = dirf( ct_r(j,2),ct_r(j,1) ); 
        p_vec(j,1) = cos( theta_temp ); 
        p_vec(j,2) = sin( theta_temp ); 
    end
    p_sign(:,1) = sign( p_asym(:,1) ); 
    p_vec = p_vec.*repmat( p_sign,[1 2] ); 
    % toc
    Plty{i,1} = [ct_1 p_vec p_asym ct{i,2}]; 
end



%% Save data
out_p.ct = ct; 
out_p.Plty = Plty; 
out_p.p_cell = p_cell; 
save( [routeOut '\Data_polarity.mat'],'out_p' )


% ========================
% end
% ========================




