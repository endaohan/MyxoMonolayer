% ***********************************************************************
% PV_finddefect
%
% Find defects using the BFld images and plot them out
% Based on bfld_finddefect.m in folder TFM_Analysis. 
% Written by Endao Han, Version 1, 8/15/2022
% 
% Input: 
% I - images to use
% dataOut - output route
% ***********************************************************************

function [dirf_all,ADefs,pdList,ndList,lays] = PV_finddefect( I,dataOut )

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


%% Setup basic parameters
dImg = [ size(I,1) size(I,2) ]; 
N_img = size( I,3 ); 

% For storing the defects
ADefs = cell(N_img,1); 
% Katie's colormap for director
% mycmap = kc_orientcmap; 

% x and y
x = 1:1:dImg(2); 
y = 1:1:dImg(1); 

% Only consider defects that last longer than this (frames)
df_length = 4; 
% Lost frames allowed
df_lost = 0; 

% Find the layer number: lays = 0 are holes
lays = find_layers( I ); 



%% Find defect
disp('Finding defects ...')
dirf_all = nan( dImg(1),dImg(2),N_img ); 

% hd = figure('Unit','Inch','Position',[2 2 12 5]); 
for i = 1:1:N_img
    % imgName = sprintf( 'Img_%04d',i ); 
    I_his = adapthisteq( rescale( I(:,:,i) ) ); 
    I_his = imadjust( I_his ); 
    % I_his = rescale( I(:,:,i) ); 
    
    % Calculate the director field and order parameter
    dirf = kc_dfield_noSave( I_his ); 
    S = kc_nemorderfield_noSave( dirf ); 
    
    % Find the defects
    try
        adefs = kc_finddefects_noSave( dirf,S,lays(:,:,i),dImg ); 
    catch
        adefs = []; 
        adefs.x = 0; 
        adefs.y = 0; 
        adefs.q = 0; 
        adefs.d = [0 0]; 
        adefs.iso = 0; 
    end
    ADefs{i,1} = adefs; 
    %}
    
    % dirf(lays(:,:,i) == 0) = NaN; 
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
    plot( ndList{i,1}.x,ndList{i,1}.y,'g.-' )
    % text( ndList{i,1}.x(1)+20,ndList{i,1}.y(1),num2str(i),'Color','g' )
end
hold off
set(gca,'YDir','reverse')
axis equal
axis([0 dImg(2) 0 dImg(1)]); 
xlabel('X (pixel)')
ylabel('Y (pixel)')
text( 800,100,'+1/2','Color','r' )
text( 800,150,'-1/2','Color','g' )
box on


% Save defect information
if ~isempty(dataOut)
    print( hd_defect,[dataOut '\DefectTracking.png'],'-dpng','-r200' )
    save( [dataOut '\dList.mat'],'pdList','ndList','lays' );
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



