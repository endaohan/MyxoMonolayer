%*************************************************************************
% Disp = TFM_disp_new(route, Parameter, flag_med)
% Based on TFM_disp_hl.m. 
% Written by Endao Han, V1, 2021/1/26
% Separate the high and low components of the displacement when doing the
% analysis, but only output the total displacement in the end. 
% Input: 
% route - input route, down to which channel to use, e.g. Beads_1
% Parameter - the Parameter structure in the TFM analysis. 
% flag_med - Use median filter (1) or not (0)
%*************************************************************************

% For debugging
%{
clear
close all
 
folder = 'D:\DATA_Confocal\20210126\PAA430Pa_Beads100nm_noCover_Coat_TimeSeries_Depth_60x_1x_1'; 

% route = [folder '\Beads_1']; 
route = [folder '\Beads_1']; 
load( [folder '\Parameters.mat'] ); 

% Use median filter on the displacement field? 1 - Yes, 0 - No. 
flag_med = 1; 
%}

function Disp = TFM_disp_new( route, Parameter, flag_med )
tic

if ~exist('flag_med','Var')
    flag_med = 0; 
end

% Parameters for making plots
S_max = 0.2*1E-6;   % maximum value in colormap
vScale = 1*1E8;     % Vector pre-factor
% Frame rate of the output videos (frames per second)
fr = Parameter.fr; 

% Set input routes
route_TFM = [route '\TFM_Data']; 
% route_Mask = [route '\Masks']; 

% Set output routes
routeOut = [route '\DispField']; 
% routeOut = [folder '\' fileList{ind_file} '\Beads_1\DispField_Ref1']; 
if ~exist(routeOut,'file')
    mkdir(routeOut)
    % mkdir([routeOut '\DispField_Data'])
    % mkdir([routeOut '\DispField_Images'])
    % mkdir([routeOut '\DispLF_Images'])
end

% List of file names
PIVList_TFM = cell( Parameter.nRef,1 ); 
N_TFM = zeros( Parameter.nRef,1 ); 
for j = 1:1:Parameter.nRef
    PIVList_TFM{j,1} = FUN_FileName( [route_TFM '\TFM_Data_Ref' ... 
        num2str(Parameter.iRef(j))], '.mat','PIV_','No' ); 
    N_TFM(j) = length( PIVList_TFM{j,1} ); 
end
% Total number of frames after PIV
N = Parameter.bd_TW(2)-Parameter.bd_TW(1); 



%% Get information from the first frame
% TFM PIV data
load([route_TFM '\TFM_Data_Ref' num2str(Parameter.iRef(1)) ... 
    '\' PIVList_TFM{1,1}{1,1} '.mat'],'Vectors'); 
% 
X = unique(Vectors.x); 
Y = unique(Vectors.y); 
% [X_mesh, Y_mesh] = meshgrid(X,Y); 
D = [length(X) length(Y)]; 

D_drift = cell(Parameter.nRef,3); 


%% Reshape data
for i = 1:1:Parameter.nRef
    % disp(i)
    % Displacement in X and Y directions in all frames
    X_tot = zeros( D(2),D(1),N_TFM(i) ); 
    Y_tot = zeros( D(2),D(1),N_TFM(i) ); 
    % Correlation coefficient 
    C_m = zeros( D(2),D(1),N_TFM(i) ); 
    
    for j = 1:1:N_TFM(i)
        load([route_TFM '\TFM_Data_Ref' num2str(Parameter.iRef(i)) ... 
            '\' PIVList_TFM{i}{j} '.mat'],'Vectors'); 
        % data for displacement field
        % data = [Vectors.x.' Vectors.y.' Vectors.u.' Vectors.v.']; 
        corrC = Vectors.correlationCoefficient; 
        
        % Reshape the velocity data into figure-sized matrices
        U_m = reshape( Vectors.u.',D ).';
        V_m = reshape( Vectors.v.',D ).';
        % S_m = sqrt(U_m.^2+V_m.^2); 
        C_m(:,:,j) = reshape( corrC,D ).'; 

        % Copy data to the total
        X_tot(:,:,j) = U_m; 
        Y_tot(:,:,j) = V_m; 
    end
        
    % Store X and Y displacements in D_drift
    D_drift{i,1} = X_tot; 
    D_drift{i,2} = Y_tot; 
    D_drift{i,3} = rescale(C_m); 
    
end


%% Shift to overlap different sections
% ind_0Hz = Parameter.bd_TW(1):1:Parameter.bd_TW(2)-1; 
% x_drift = zeros( D(2),D(1),length(ind_0Hz) ); 
% y_drift = zeros( D(2),D(1),length(ind_0Hz) ); 

% Take sections from each reference frame and stitch them together
for i = 1:1:Parameter.nRef
    if i == 1
        D_drift{i,1} = D_drift{i,1}; 
        D_drift{i,2} = D_drift{i,2}; 
    else
        % Indices of the overlapped images in D_drift
        i1 = ( Parameter.tRange(i,2):1:Parameter.tRange(i-1,3)-1 )-Parameter.tRange(i-1,2)+1;
        i2 = ( Parameter.tRange(i,2):1:Parameter.tRange(i-1,3)-1 )-Parameter.tRange(i,2)+1; 
        % Mean of the overlapped region
        x1 = trimmean( D_drift{i-1,1}(:,:,i1),10,3 ); 
        x2 = trimmean( D_drift{i,1}(:,:,i2),10,3 ); 
        y1 = trimmean( D_drift{i-1,2}(:,:,i1),10,3 ); 
        y2 = trimmean( D_drift{i,2}(:,:,i2),10,3 ); 
        % Get rid of the mean
        D_drift{i,1} = D_drift{i,1}-x2+x1; 
        D_drift{i,2} = D_drift{i,2}-y2+y1; 
    end
end



%% Calculate the mean displacement
U_mean = zeros( D(2),D(1),N ); 
V_mean = zeros( D(2),D(1),N ); 
U_std = zeros( D(2),D(1),N ); 
V_std = zeros( D(2),D(1),N ); 
NFrm = zeros(N,1); 
C_mean = zeros( D(2),D(1),N ); 
% i_Uex = zeros( D(2),D(1),N_tot ); 
% i_Vex = zeros( D(2),D(1),N_tot ); 

for j = 1:1:N
    frmName = sprintf('PIV_VectorData%06d',j);

    % Find the number of reference frames this frame is related to
    % 1 - reference frame index; 2 - index of this frame in PIVList_TFM
    ind_frm = []; 
    for k = 1:1:Parameter.nRef
        indStr = strcmp(PIVList_TFM{k,1},frmName); 
        if sum( indStr ) == 1
            ind_frm = [ind_frm; k find(indStr == 1)]; 
        end
    end
    nFrm = size(ind_frm,1); 
    NFrm(j) = nFrm; 
    
    % Calculate the mean displacement in each frame
    U_m = zeros( D(2),D(1),nFrm ); 
    V_m = zeros( D(2),D(1),nFrm ); 
    C_m = zeros( D(2),D(1),nFrm ); 
    % load data in the range to be analyzed
    for k = 1:1:nFrm
        U_m(:,:,k) = D_drift{ind_frm(k,1),1}(:,:,ind_frm(k,2)); 
        V_m(:,:,k) = D_drift{ind_frm(k,1),2}(:,:,ind_frm(k,2)); 
        C_m(:,:,k) = D_drift{ind_frm(k,1),3}(:,:,ind_frm(k,2)); 
    end
    % Mean correlation coefficient
    C_mean(:,:,j) = mean(C_m,3); 
    
    % Calculate the mean while removing outliers
    iOut = isoutlier(U_m,'median',3); 
    U_m( iOut==1 ) = NaN; 
    iOut = isoutlier(V_m,'median',3); 
    V_m( iOut==1 ) = NaN; 
    
    U_mean(:,:,j) = mean( U_m,3,'omitnan' ); 
    V_mean(:,:,j) = mean( V_m,3,'omitnan' ); 
    U_std(:,:,j) = std( U_m,0,3,'omitnan' ); 
    V_std(:,:,j) = std( V_m,0,3,'omitnan' ); 
    
end

% Make the first frame the reference frame
U_mean = U_mean-repmat( U_mean(:,:,1),1,1,N ); 
V_mean = V_mean-repmat( V_mean(:,:,1),1,1,N ); 

% Median filter
if flag_med == 1
    U_mean = medfilt3( U_mean,[3 3 1] ); 
    V_mean = medfilt3( V_mean,[3 3 1] ); 
end




%% Remove noisy data
% Use three parameters to decide which regions are to be excluded in the
% final force map: noise to signal ratio, magnitude of displacement, and 
% the correlation coefficient from PIV C_mean. 
Dm = sqrt(U_mean.^2+V_mean.^2); 
% NtS = sqrt(U_std.^2+V_std.^2)./Dm; 
NtS = sqrt(U_std.^2+V_std.^2); 
D_med = median( Dm,'all' ); 
NtS(isinf(NtS)) = NaN; 
NtS_thresh = mean(NtS,'all','omitnan')+std(NtS,0,'all','omitnan'); 
C_thresh = mean(C_mean,'all','omitnan')-3*std(C_mean,0,'all','omitnan'); 
ind_ex = zeros( size(U_mean) ); 

% Excluded pixels are labelled as 1
ind_ex( (NtS>NtS_thresh)&(Dm>D_med)|(C_mean<C_thresh) ) = 1; 
ind_ex = medfilt3( ind_ex,[5 5 5] ); 
ind_ex = logical(ind_ex); 

% U_mean( ind_ex == 1 ) = NaN; 
% V_mean( ind_ex == 1 ) = NaN; 
% U_lf( ind_ex == 1 ) = 0; 
% V_lf( ind_ex == 1 ) = 0; 

% U_mean = fillmissing( U_mean,'movmean',3 );  
% V_mean = fillmissing( V_mean,'movmean',3 );  
% U_lf = fillmissing( U_lf,'movmedian',5 );  
% V_lf = fillmissing( V_lf,'movmedian',5 );  





%% Slow drift
% Use gaussian filter to separate frequencies
win_s = 0.5; 
win_t = 50;  
if Parameter.N_img < 40
    win_end = 5; 
else
    win_end = 20; 
end

%{
% U_lf = medfilt3( U_mean,[3,3,3] ); 
% V_lf = medfilt3( V_mean,[3,3,3] ); 
U_lf = imgaussfilt3(U_mean,[win_s win_s win_t],'FilterSize',[1 1 2*ceil(2*win_t)+1]);
V_lf = imgaussfilt3(V_mean,[win_s win_s win_t],'FilterSize',[1 1 2*ceil(2*win_t)+1]);
%}

% Pad ends and gaussian filter
% Displacement in x
U_first = mean(U_mean(:,:,1:win_end),3); 
U_end = mean(U_mean(:,:,end-win_end+1:end),3); 
U_mean = cat(3,U_first,U_mean,U_end); 
U_lf = imgaussfilt3(U_mean,[win_s win_s win_t],'FilterSize',[1 1 2*ceil(2*win_t)+1]);
U_lf(:,:,[1 end]) = []; 
U_mean(:,:,[1 end]) = []; 

% Displacement in y
V_first = mean(V_mean(:,:,1:win_end),3); 
V_end = mean(V_mean(:,:,end-win_end+1:end),3); 
V_mean = cat(3,V_first,V_mean,V_end); 
V_lf = imgaussfilt3(V_mean,[win_s win_s win_t],'FilterSize',[1 1 2*ceil(2*win_t)+1]);
V_lf(:,:,[1 end]) = []; 
V_mean(:,:,[1 end]) = []; 

U_hf = U_mean-U_lf; 
V_hf = V_mean-V_lf; 

% Use high pass filter
%{
U_hf = zeros( D(2),D(1),N ); 
V_hf = zeros( D(2),D(1),N ); 

for j = 1:1:D(1)
    U_temp = reshape( U_mean(:,j,:),[D(2) N] ); 
    V_temp = reshape( V_mean(:,j,:),[D(2) N] ); 
    U_hf(:,j,:) = highpass( U_temp.',1/1E3,1/Parameter.tscl,'Steepness',0.95 ).'; 
    V_hf(:,j,:) = highpass( V_temp.',1/1E3,1/Parameter.tscl,'Steepness',0.95 ).'; 
end
    
U_lf = U_mean-U_hf; 
V_lf = V_mean-V_hf; 
%}


%{
figure; 
dif_p = floor(N_tot/5); 
ind_p = 1:dif_p:N_tot; 
for i = 1:1:length(ind_p)
    subplot(2,length(ind_p),i)
    imagesc(U_lf(:,:,i))
    caxis([0 S_max])
    subplot(2,length(ind_p),length(ind_p)+i)
    imagesc(V_lf(:,:,i))
    caxis([0 S_max/10])
end
colorbar
%}


% Plot stdev of displacement before and after removing the low frequency
% mode. 
pstdev.U_orig = reshape( std(U_mean,0,[1 2]),[N 1] )*1E6; 
pstdev.U_hf = reshape( std(U_hf,0,[1 2]),[N 1] )*1E6; 
pstdev.V_orig = reshape( std(V_mean,0,[1 2]),[N 1] )*1E6; 
pstdev.V_hf = reshape( std(V_hf,0,[1 2]),[N 1] )*1E6; 

hd_stdev = figure; 
hold on
plot( pstdev.U_orig,'b--' )
plot( pstdev.U_hf,'b-' )
plot( pstdev.V_orig,'m--' )
plot( pstdev.V_hf,'m-' )
% plot([21 21],[0 0.1],'k--')
% plot([length(a)-20 length(a)-20],[0 0.1],'k--')
hold off
xlabel('Frame')
ylabel('Stdev displacement (\mum)')
legend('x original','x high frequency','y original','y high frequency', ... 
    'Location','northwest')


%{
figure; 
subplot(1,2,1)
imagesc( max(U_lf,[],3)-min(U_lf,[],3) ); 
colorbar
subplot(1,2,2)
imagesc( max(V_lf,[],3)-min(V_lf,[],3) ); 
colorbar
%}


%% FFT
[f,~] = FUN_FFT_v2(1/15,reshape(U_hf(1,1,:),[N 1]),0); 
Nf = length(f); 
P1 = zeros( D(2),D(1),Nf ); 
P2 = zeros( D(2),D(1),Nf ); 
P3 = zeros( D(2),D(1),Nf ); 

for i = 1:1:D(2)
    for j = 1:1:D(1)
        [~,P1(i,j,:)] = FUN_FFT_v2(1/15,reshape(U_hf(i,j,:),[N 1]),0); 
        [~,P2(i,j,:)] = FUN_FFT_v2(1/15,reshape(U_lf(i,j,:),[N 1]),0); 
        [~,P3(i,j,:)] = FUN_FFT_v2(1/15,reshape(U_mean(i,j,:),[N 1]),0); 
    end
end

pfft.f = f; 
pfft.P_hf = reshape( mean(P1,[1 2]),[length(f) 1] ); 
pfft.P_lf = reshape( mean(P2,[1 2]),[length(f) 1] ); 
pfft.P = reshape( mean(P3,[1 2]),[length(f) 1] ); 

hd_fft = figure; 
hold on
% plot( pfft.f,10*log10(pfft.P_hf),'-' )
% plot( pfft.f,10*log10(pfft.P_lf),'-' )
% plot( pfft.f,10*log10(pfft.P),'-' )
plot( pfft.f,pfft.P_hf,'-' )
plot( pfft.f,pfft.P_lf,'-' )
plot( pfft.f,pfft.P,'-' )
hold off
ylabel('Amplitude spectrum')
% ylabel('PSD (dB/Hz)')
xlabel('Frequency (Hz)')
legend('High frequency','Low frequency','Total', ... 
    'Location','northeast')




%% Calculate the zero frequency mode
%{
D_0Hz = cell(1,2); 
% D_0Hz{1,1} = mean( U_hf,3,'omitnan' ); 
% D_0Hz{1,2} = mean( V_hf,3,'omitnan' ); 
D_0Hz{1,1} = 0; 
D_0Hz{1,2} = 0; 

U_0Hz = U_hf-D_0Hz{1,1}; 
V_0Hz = V_hf-D_0Hz{1,2}; 
    
save([routeOut '\Drift0Hz.mat'],'D_0Hz'); 
%} 

% Make a plot to show that the shifting is good
hd_shift = figure('Unit','Inch','Position',[2 2 6 8]); 
% ip = [60 45]; 
ip = round( [D(2) D(1)]/2 ); 
subplot(3,1,1)
hold on
for i = 1:1:size(D_drift,1)
    t = Parameter.tRange(i,2):1:Parameter.tRange(i,3)-1; 
    f = reshape( D_drift{i,1}(ip(1),ip(2),:),[length(t) 1] )*1E6; 
    plot( t,f,'b-' )
    f = reshape( D_drift{i,2}(ip(1),ip(2),:),[length(t) 1] )*1E6; 
    plot( t,f,'r-' )
end
plot( reshape(U_mean(ip(1),ip(2),:)*1E6,[N 1]),'k-' )
plot( reshape(V_mean(ip(1),ip(2),:)*1E6,[N 1]),'k-' )
hold off
box on
ylim([-0.1 0.1])
xlabel('Frame')
ylabel('Displacement (\mum)')
title('Original displacement')

subplot(3,1,2)
hold on
plot( reshape(U_lf(ip(1),ip(2),:)*1E6,[N 1]),'b-' )
plot( reshape(V_lf(ip(1),ip(2),:)*1E6,[N 1]),'r-' )
hold off
box on
ylim([-0.1 0.1])
xlabel('Frame')
ylabel('Displacement (\mum)')
title('Low frequency mode')
legend('x','y','Location','northwest')

subplot(3,1,3)
hold on
plot( reshape(U_hf(ip(1),ip(2),:)*1E6,[N 1]),'b-' )
plot( reshape(V_hf(ip(1),ip(2),:)*1E6,[N 1]),'r-' )
hold off
box on
ylim([-0.05 0.05])
xlabel('Frame')
ylabel('Displacement (\mum)')
title('High frequency mode')




%% Output data
% Only include the origianl displacement
% iOut = [21 40 15 40]; 
iOut = [1 D(2) 1 D(1)]; 
Disp.X = X(iOut(3):iOut(4)); 
Disp.Y = Y(iOut(1):iOut(2)); 
Disp.U = U_mean( iOut(1):iOut(2), iOut(3):iOut(4),: ); 
Disp.V = V_mean( iOut(1):iOut(2), iOut(3):iOut(4),: ); 
Disp.ind_ex = ind_ex( iOut(1):iOut(2), iOut(3):iOut(4),: );  




%% Plot correlation time and length scales
% Spatial correlation
Cs = corr_spatial_vec( U_hf,V_hf ); 
Cs_lf = corr_spatial_vec( U_lf,V_lf ); 
db = bin_ave( Cs.R,Cs.C,100 ); 
db = [0 0 Cs.C(Cs.R==0) 0 1; db]; 
n_grid = X(2)-X(1); 
% Temporal correlation
Ct = corr_time_vec( U_hf,V_hf ); 
Ct_lf = corr_time_vec( U_lf,V_lf ); 


% Make plot
hd_corr = figure; 
subplot(1,2,1)
plot( Cs_lf.R*n_grid*Parameter.lscl(1),Cs_lf.C,'r.' )
hold on
plot( Cs.R*n_grid*Parameter.lscl(1),Cs.C,'b.' )
% plot( db(:,2)*n_grid*Parameter.lscl(1),db(:,3),'k-' )
% plot( xq*n_grid*Parameter.lscl(1),yq,'r-' )
% plot( xplot,pf(1)*xplot.^pf(2).*exp( pf(3)*xplot ),'m-' )
% plot( xplot,F(pf,xplot),'r-','LineWidth',2 )
% plot( xplot,g(pg,xplot),'g-' )
hold off
xlim([0 50])
ylim([-0.2 1])
xlabel('r (\mum)')
ylabel('Autocorrelation')
legend('Low frequency','High frequency', ... 
    'Location','northeast')

subplot(1,2,2)
plot( Ct.T*Parameter.tscl,Ct.C,'bo' )
hold on
plot( Ct_lf.T*Parameter.tscl,Ct_lf.C,'ro' )
plot( [0 max(Ct.T*Parameter.tscl)],[0 0],'k--' )
hold off
ylim([-0.2 1])
xlabel('Time (s)')
ylabel('Autocorrelation')

pcorr.r = Cs.R; 
pcorr.Cr = Cs.C; 
pcorr.lscl = n_grid*Parameter.lscl(1); 
pcorr.t = Ct.T; 
pcorr.tscl = Parameter.tscl; 
pcorr.Ct = Ct.C; 
pcorr.n = Ct.n; 




%% Save videos for high and low frequency modes 

% Displacement of total displacement
S_m = sqrt( U_mean.^2 + V_mean.^2 ); 

videoName = 'Disp_total'; 
v = VideoWriter([routeOut '\' videoName '.mp4'],'MPEG-4'); 
% v = VideoWriter([routeOut '\' videoName '.avi'],'Uncompressed AVI'); 
v.FrameRate = fr; 
disp('Making displacement video...')

hd = figure('units','inch','position',[1,0,5,6]); 
% Start making video
open(v)
% Loop frame by frame
for j = 1:1:N
    imagesc( X,Y,S_m(:,:,j) ); 
    hold on
    quiver(X(1:2:end),Y(1:2:end), ... 
        U_mean(1:2:end,1:2:end,j)*vScale, ... 
        -V_mean(1:2:end,1:2:end,j)*vScale, ... 
        'k','AutoScale','off')
    % Label the points that have poor data quality
    [yn,xn] = find( ind_ex(:,:,j)==1 ); 
    plot( X(xn),Y(yn),'ro','LineWidth',0.6)
    hold off
    
    ax = gca;
    ax.YDir = 'reverse'; 
    axis equal
    axis([0 1000 0 1400])
        
    caxis([0 S_max])
    colormap('cool')
    colorbar
    
    ImgVideo = getframe(hd); 
    writeVideo( v,ImgVideo ); 
    
    clf
end
close(v)


% Displacement of high frequency mode
S_m = sqrt( U_hf.^2 + V_hf.^2 ); 

videoName = 'Disp_HF'; 
v = VideoWriter([routeOut '\' videoName '.mp4'],'MPEG-4'); 
% v = VideoWriter([routeOut '\' videoName '.avi'],'Uncompressed AVI'); 
v.FrameRate = fr; 
disp('Making displacement video...')

hd = figure('units','inch','position',[1,0,5,6]); 
% hd = figure('units','pixels','position',[0 0 864 1152]); 
% Start making video
open(v)

% Loop frame by frame
for j = 1:1:N
    imagesc( X,Y,S_m(:,:,j) ); 
    hold on
    quiver(X(1:2:end),Y(1:2:end), ... 
        U_hf(1:2:end,1:2:end,j)*vScale*2, ... 
        -V_hf(1:2:end,1:2:end,j)*vScale*2, ... 
        'k','AutoScale','off')
    % Label the points that have poor data quality
    [yn,xn] = find( ind_ex(:,:,j)==1 ); 
    plot( X(xn),Y(yn),'ro','LineWidth',0.6)
    hold off
    
    ax = gca;
    ax.YDir = 'reverse'; 
    axis equal
    axis([0 1000 0 1400])
        
    caxis([0 S_max/2])
    colormap('cool')
    colorbar
    
    ImgVideo = getframe(hd); 
    writeVideo(v, ImgVideo); 
    
    clf
end
close(v)


% Displacement of low frequency mode
S_m = sqrt( U_lf.^2 + V_lf.^2 ); 

videoName = 'Disp_LF'; 
% v = VideoWriter([routeOut '\' videoName '.avi'],'Uncompressed AVI'); 
v = VideoWriter([routeOut '\' videoName '.mp4'],'MPEG-4'); 
v.FrameRate = fr; 
disp('Making displacement video...')

hd_lf = figure('units','inch','position',[1,0,5,6]); 
% hd_lf = figure('units','pixels','position',[0 0 690 920]); 
% Start making video
open(v)

% Loop frame by frame
for j = 1:1:N
    % frmName = sprintf('PIV_VectorData%06d',j);

    imagesc( X,Y,S_m(:,:,j) ); 
    hold on
    quiver(X(1:2:end),Y(1:2:end), ... 
        U_lf(1:2:end,1:2:end,j)*vScale, ... 
        -V_lf(1:2:end,1:2:end,j)*vScale, ... 
        'k','AutoScale','off')
    hold off
    
    ax = gca;
    ax.YDir = 'reverse'; 
    axis equal
    axis([0 1000 0 1400])
        
    caxis([0 S_max])
    colormap('cool')
    colorbar
    
    % saveas(hd, [routeOut '\DispField_Images\' frmName, '.tif']);
    % rint(hd_lf, [routeOut '\DispLF_Images\' frmName, '.tif'], '-dtiff', '-r200')
    % a = getframe(hd_lf); 
    % imwrite(a.cdata,[routeOut '\DispLF_Images\' frmName, '.png']); 
    ImgVideo = getframe(hd_lf); 
    writeVideo(v, ImgVideo); 
    
    clf       
end
close(v)



%% Save results
% saveas(hd, [routeOut '\DispField_Images\' frmName, '.tif']);
print(hd_shift, [routeOut '\MeanDisplacement.png'], '-dpng', '-r200')
print(hd_stdev, [routeOut '\Stdev_frames.png'], '-dpng', '-r200')
print(hd_fft, [routeOut '\FFT.png'], '-dpng', '-r200')
print(hd_corr, [routeOut '\Correlation.png'], '-dpng', '-r200')

% Output data for future reference
dstdev.U_std = U_std; 
dstdev.V_std = V_std; 
save([routeOut '\outMat.mat'], 'dstdev','pfft','pstdev','pcorr');

% Save the displacement data for force reconstruction
save([routeOut '\Disp.mat'], 'Disp');

toc
%==============================
end
%==============================


%% Local functions





