% ***********************************************************************
% BFld_finddefect_test_v7( folder,fileName,flag_plot )
%
% Find defects using the BFld images and plot them out
% Based on BFld_finddefect_test_v9.m in folder 20210128. 
% Written by Endao Han, Version 2, 8/3/2021
% 
% Input: 
% folder - folder the data in (to date)
% fileName - name of the file
% flag_plot - plot the defects (1) or not (0)
% ***********************************************************************

function [pdList,ndList] = BFld_finddefect( folder,fileName,flag_plot )

% For debugging
%{
% date = '20210126'; 
% date = '20201014'; 
% fileName = 'PAA430Pa_Beads100nm_noCover_Coat_TimeSeries_Depth_60x_1x_1'; 
% date = '20210205'; 
% fileName = 'PAA430Pa_Beads100nm_noCover_Coat_TimeSeries_Depth_60x_1x_4'; 
date = '20201006'; 
fileName = 'PAA430Pa_Beads40100nm_noCover_Coat_TimeSeries_60x_1x_4'; 
folder = ['D:\DATA_Confocal\' date]; 
%}

% Flag_plot
if ~exist('flag_plot','var')
    flag_plot = 0; 
end


%% Setup basic parameters
routeIn = [ folder '\' fileName ]; 
fpath = [routeIn '\BFld']; 
route = [routeIn '\BFld\Images']; 
% routeOut = [routeIn '\DField']; 

load([routeIn '\Parameters.mat'],'Parameter'); 
N_img = Parameter.bd_TW(2)-1; 
% N_img = 150; 

% For storing the defects
ADefs = cell(N_img,1); 
% Katie's colormap for director
% mycmap = kc_orientcmap; 

% x and y
x = 1:1:Parameter.dImg(2); 
y = 1:1:Parameter.dImg(1); 

% Only consider defects that last longer than this (frames)
df_length = 8; 
% Lost frames allowed
df_lost = 0; 



%% Read in images
I = zeros( Parameter.dImg(1),Parameter.dImg(2),N_img ); 
% I_mask = zeros( Parameter.dImg(1),Parameter.dImg(2),N_img ); 

for i = 1:1:N_img
    imgName = sprintf( 'Img_%04d',i ); 
    % I(:,:,i) = double( imread([route '\' imgName '.tif']) ); 
    I(:,:,i) = imread([route '\' imgName '.tif']);  
end
% Find the layer number: lays = 0 are holes
lays = find_layers( I ); 



%% Find defect
disp('Finding defects ...')

dirf_all = zeros( Parameter.dImg(1),Parameter.dImg(2),N_img ); 

% hd = figure('Unit','Inch','Position',[2 2 12 5]); 
for i = 1:1:N_img
    % imgName = sprintf( 'Img_%04d',i ); 
    % I(:,:,i) = imread( [route '\' imgName '.tif'] ); 
    I_temp = rescale( I(:,:,i) );
    I_his = adapthisteq( I_temp ); 
    I_his = imadjust(I_his); 
    
    % Calculate the director field and order parameter
    dirf = kc_dfield( fpath,I_his,i ); 
    kc_nemorderfield( fpath,i,Parameter.dImg ); 
    
    % Find the defects
    try
        adefs = kc_finddefects( fpath,i,lays(:,:,i),Parameter.dImg ); 
    catch
        adefs = []; 
        adefs.x = 0; 
        adefs.y = 0; 
        adefs.q = 0; 
        adefs.d = [0 0]; 
        adefs.iso = 0; 
    end
    ADefs{i,1} = adefs; 
    
    dirf(lays(:,:,i) == 0) = NaN; 
    dirf_all(:,:,i) = dirf; 
    % S = S.*double(lays(:,:,i)); 
end



%% Separate positive and negative defects
% Pos and neg defects frame by frame
Apos = cell(N_img,1); 
Aneg = cell(N_img,1); 

disp('Separate defects...')
for i = 1:1:N_img
    ind_pn = nan( length(ADefs{i,1}),1 );
    for j = 1:1:length(ADefs{i,1})        
        if ADefs{i,1}(j).q > 0.1
            ind_pn(j) = 1; 
        elseif ADefs{i,1}(j).q < -0.1
            ind_pn(j) = 0; 
        else
            ind_pn(j) = NaN; 
        end
    end
    
    def_temp = ADefs{i,1}; 
    % def_temp(ind_pn == 0) = []; 
    Apos{i,1} = def_temp( ind_pn == 1 ); 
    
    def_temp = ADefs{i,1}; 
    % def_temp(ind_pn == 1) = []; 
    Aneg{i,1} = def_temp( ind_pn == 0 ); 
end



%% Defect tracking
% Track defects using Apos and Aneg. 
% The output contain two columns: 1 - data structure; 2 - starting and
% ending frame number. 
% In the structure, id - id number for the defect; x,y - x,y positions of
% the defect in pixel; i - frame. 
disp('Tracking defects')
pdList = track_defect( Apos,df_lost,df_length ); 
ndList = track_defect( Aneg,df_lost,df_length ); 


%% Assign director d to every defect
pdList = assign_director( pdList,Apos ); 
ndList = assign_director( ndList,Aneg ); 


        
%% Make plot of all the defects
% Locations of all the positive and negative defects identified
hd_defect = figure; 
hold on
for i = 1:1:length(pdList)
    plot( pdList{i,1}.x,pdList{i,1}.y,'r.-' )
    % text( pdList{i,1}.x(1)+20,pdList{i,1}.y(1),num2str(i),'Color','r' )
end
for i = 1:1:length(ndList)
    plot(ndList{i,1}.x,ndList{i,1}.y,'g.-')
    % text( ndList{i,1}.x(1)+20,ndList{i,1}.y(1),num2str(i),'Color','g' )
end
hold off
set(gca,'YDir','reverse')
axis equal
axis([0 Parameter.dImg(2) 0 Parameter.dImg(1)]); 
xlabel('X (pixel)')
ylabel('Y (pixel)')
text( 800,100,'+1/2','Color','r' )
text( 800,150,'-1/2','Color','g' )
box on


% Save defect information
dataOut = [fpath '\analysis']; 
if ~exist( dataOut,'file' )
    mkdir( dataOut ); 
end
print( hd_defect,[dataOut '\DefectTracking.png'],'-dpng','-r200' )
save( [dataOut '\dList.mat'],'pdList','ndList','lays' );



%% Plot the tracking on BFld and stress map and make videos
if flag_plot
    load([folder '\' fileName '\Force\Trac.mat'],'Trac'); 
    F_sm = sqrt( Trac.Fx.^2+Trac.Fy.^2 ); 
    % F_sm = sqrt( Trac.Fx_hf.^2+Trac.Fy_hf.^2 ); 

    % Length scale for the defect label
    vlscl = 40; 
    % Length scale for the traction vectors
    vScale = 0.5; 
    % rotation matrix
    Rm = [cos(2*pi/3) -sin(2*pi/3); sin(2*pi/3) cos(2*pi/3)]; 

    % Output route
    routeP = [fpath '\analysis\Stress']; 
    if ~exist( routeP,'file' )
        mkdir( routeP ); 
    end
    
    % Make the figures
    hd = figure('Unit','Inch','Position',[2 2 10 6]); 
    for i = 1:1:N_img
        % ax(1) = subplot(1,2,1,'align'); 
        ax(1) = subplot('Position',[0.05 0.07 0.4 0.9]); 
        imagesc( I(:,:,i) ); 
        colormap(ax(1),'gray')
        hold on 
        imagesc( x,y, ones(Parameter.dImg)*0.5, ... 
            'AlphaData',(1-lays(:,:,i))*0.3 )
        for j = 1:1:length(pdList)
            if ismember( i,pdList{j,1}.i )
                ip = find( pdList{j,1}.i == i ); 
                plot( pdList{j,1}.x(ip),pdList{j,1}.y(ip), ... 
                    'ro','MarkerSize',4, ... 
                    'MarkerFaceColor','r' )
                r = pdList{j,1}.d(ip,:).'*vlscl; 
                plot( [pdList{j,1}.x(ip) pdList{j,1}.x(ip)+r(1)], ... 
                    [pdList{j,1}.y(ip) pdList{j,1}.y(ip)+r(2)], ... 
                    'r-','LineWidth',1.5 )
                text( pdList{j,1}.x(ip)+15,pdList{j,1}.y(ip)-10, ... 
                    num2str(pdList{j,1}.id), ... 
                    'FontSize',12,'Color',[1 0.4 0.2] )
            end
        end
    
        for j = 1:1:length(ndList)
            if ismember(i,ndList{j,1}.i)
                ip = find( ndList{j,1}.i == i ); 
                plot( ndList{j,1}.x(ip),ndList{j,1}.y(ip), ... 
                    'o','MarkerSize',4, ... 
                    'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1] )
                r = ndList{j,1}.d(ip,:).'*vlscl; 
                plot([ndList{j,1}.x(ip) ndList{j,1}.x(ip)+r(1)], ... 
                    [ndList{j,1}.y(ip) ndList{j,1}.y(ip)+r(2)], ... 
                    '-','LineWidth',1.5,'Color',[0 0 1])
                r = Rm*r; 
                plot([ndList{j,1}.x(ip) ndList{j,1}.x(ip)+r(1)], ... 
                    [ndList{j,1}.y(ip) ndList{j,1}.y(ip)+r(2)], ... 
                    '-','LineWidth',1.5,'Color',[0 0 1])
                r = Rm*r; 
                plot([ndList{j,1}.x(ip) ndList{j,1}.x(ip)+r(1)], ... 
                    [ndList{j,1}.y(ip) ndList{j,1}.y(ip)+r(2)], ... 
                    '-','LineWidth',1.5,'Color',[0 0 1])
                text( ndList{j,1}.x(ip)+15,ndList{j,1}.y(ip)-10, ... 
                    num2str(ndList{j,1}.id), ... 
                    'FontSize',12,'Color',[0.2 0.4 1] )

            end
        end
        hold off
        axis equal
        axis( [0 Parameter.dImg(2) 0 Parameter.dImg(1)] ); 
        xlabel('x (pixel)')
        ylabel('y (pixel)')

        % ax(2) = subplot(1,2,2,'align'); 
        ax(2) = subplot( 'Position',[0.5 0.07 0.5 0.9] ); 
        imagesc( Trac.X,Trac.Y,F_sm(:,:,i) )
        colormap(ax(2),'summer')
        caxis([0 50])
        hold on 

        imagesc( x,y, ones(Parameter.dImg)*0.5, ... 
            'AlphaData',(1-lays(:,:,i))*0.3 )
        % imagesc( x,y, zeros(Parameter.dImg), ... 
        %     'AlphaData',1-S_all(:,:,i) )

        for j = 1:1:length(pdList)
            if ismember( i,pdList{j,1}.i )
                ip = find( pdList{j,1}.i == i ); 
                plot(pdList{j,1}.x(ip),pdList{j,1}.y(ip), ... 
                    'ro','MarkerSize',4, ... 
                    'MarkerFaceColor','r')
                r = pdList{j,1}.d(ip,:).'*vlscl; 
                plot([pdList{j,1}.x(ip) pdList{j,1}.x(ip)+r(1)], ... 
                    [pdList{j,1}.y(ip) pdList{j,1}.y(ip)+r(2)], ... 
                    'r-','LineWidth',1.5)
            end
        end
    
        for j = 1:1:length(ndList)
            if ismember(i,ndList{j,1}.i)
                ip = find( ndList{j,1}.i == i ); 
                plot(ndList{j,1}.x(ip),ndList{j,1}.y(ip), ... 
                    'o','MarkerSize',4, ... 
                    'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1])
                r = ndList{j,1}.d(ip,:).'*vlscl; 
                plot([ndList{j,1}.x(ip) ndList{j,1}.x(ip)+r(1)], ... 
                    [ndList{j,1}.y(ip) ndList{j,1}.y(ip)+r(2)], ... 
                    '-','LineWidth',1.5,'Color',[0 0 1])
                r = Rm*r; 
                plot([ndList{j,1}.x(ip) ndList{j,1}.x(ip)+r(1)], ... 
                    [ndList{j,1}.y(ip) ndList{j,1}.y(ip)+r(2)], ... 
                    '-','LineWidth',1.5,'Color',[0 0 1])
                r = Rm*r; 
                plot([ndList{j,1}.x(ip) ndList{j,1}.x(ip)+r(1)], ... 
                    [ndList{j,1}.y(ip) ndList{j,1}.y(ip)+r(2)], ... 
                    '-','LineWidth',1.5,'Color',[0 0 1])

            end
        end
        
        % Regions pushed down
        [yd,xd] = find( isnan(F_sm(:,:,i)) ); 
        plot( Trac.X(xd),Trac.Y(yd),'o','Color',[1 0.8 0],'LineWidth',0.5 )
        
        
        % Coordinates of dirf is defined the same as image, positive y is down.
        quiver( x(1:10:end), y(1:10:end), ... 
            cos(dirf_all(1:10:end,1:10:end,i))*vScale*20, ... 
            sin(dirf_all(1:10:end,1:10:end,i))*vScale*20, ... 
            'Color','w', ... 
            'ShowArrowHead','off', ... 
            'AutoScale','off', 'LineWidth',0.2)
        % Coordinates of Trac is defined differently, positive y is up.
        quiver( Trac.X(1:1:end,1:1:end), ... 
            Trac.Y(1:1:end,1:1:end), ... 
            Trac.Fx(1:1:end,1:1:end,i)*vScale, ... 
            -Trac.Fy(1:1:end,1:1:end,i)*vScale, ... 
            'Color','k', ... 
            'AutoScale','off', 'LineWidth',0.6)
        %{
        quiver( Trac.X(1:2:end,1:2:end), ... 
            Trac.Y(1:2:end,1:2:end), ... 
            Trac.Fx_hf(1:2:end,1:2:end,i)*vScale, ... 
            -Trac.Fy_hf(1:2:end,1:2:end,i)*vScale, ... 
            'Color','k', ... 
            'AutoScale','off', 'LineWidth',0.6)
        %}
        hold off
        axis equal
        axis([0 Parameter.dImg(2) 0 Parameter.dImg(1)]); 
        colorbar
        xlabel('x (pixel)')
        ylabel('y (pixel)')

        imgName = sprintf('Img_%04d',i); 
        % saveas(hd,[routeP '\' imgName '.tif']); 
        print(hd, [routeP '\' imgName '.png'],'-dpng','-r200')
        clf
        % close all
    end



    % Make video
    videoOut = [fpath '\analysis']; 
    % BFld images
    videoIn = [fpath '\analysis\Stress']; 
    videoName = 'Defects_Stress'; 
    FUN_ImgToVideo( videoIn, 'png', videoOut,videoName,7,'avi' ); 

end



% ========
end
% ========




%% Local functions

% Track the defects
function dtr = track_defect( Adef,df_lost,df_length )
    % disp('Tracking defects')
    
    % Maximum displacement (pixels)
    disp_max = 40; 
    
    % Parameter for track
    param.mem = df_lost;
    param.good = df_length;
    param.dim = 2;
    param.quiet = 0;

    pp = []; 
    for j = 1:1:length(Adef)
        pp_temp = zeros( length(Adef{j,1}),3 ); 
        pp_temp(:,1) = vertcat( Adef{j,1}.x ); 
        pp_temp(:,2) = vertcat( Adef{j,1}.y ); 
        pp_temp(:,3) = j; 
        pp = vertcat(pp,pp_temp);  
    end

    % Track using Dufresne & Blair's PTV code
    tr = track( pp,disp_max,param ); 

    % Number of defects
    n_id = max( tr(:,4) ); 
    dtr = cell(n_id,2); 
    for j = 1:1:n_id
        i_temp = find( tr(:,4) == j );
        dtr{j,1}.id = j; 
        dtr{j,1}.x = tr( i_temp,1 ); 
        dtr{j,1}.y = tr( i_temp,2 ); 
        dtr{j,1}.i = tr( i_temp,3 ); 
        dtr{j,2} = [dtr{j,1}.i(1) dtr{j,1}.i(end)]; 
    end
end


% Find the director for each defect
% And if the defect is isolated
function dList = assign_director( dList,Adef )
    nd = size(dList,1); 
    for i = 1:1:nd
        nf = length( dList{i,1}.i ); 
        d_temp = zeros(nf,2); 
        iso_temp = zeros(nf,1); 

        for j = 1:1:nf
            frm = dList{i,1}.i(j); 
            xd = vertcat( Adef{frm,1}.x ); 
            yd = vertcat( Adef{frm,1}.y ); 

            ind =  (xd == dList{i,1}.x(j))&(yd == dList{i,1}.y(j)) ; 
            d_temp(j,:) = Adef{frm,1}(ind).d; 
            iso_temp(j) = Adef{frm,1}(ind).iso; 
        end     
        dList{i,1}.d = d_temp; 
        dList{i,1}.iso = iso_temp; 
    end
end



