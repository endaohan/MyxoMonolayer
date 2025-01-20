% Use optical flow to calculate the cell velocity
function [Vx_box,Vy_box] = cal_velocity_optFlow( I_bfld,X,Y )

% Window size of Gaussian smoothing
winGauss = 32; 
% Window for local min-max filter
win_minmax = ones(6); 

N = size( I_bfld,3 );  

I_of = nan( size(I_bfld) ); 
for i = 1:1:N
    % I_temp = imadjust( method_Katie(I_proc(:,:,i),winGauss) ); 
    I_temp = imadjust( method_Katie(I_bfld(:,:,i),winGauss) ); 
    I_of(:,:,i) = FUN_localMinMax( I_temp.^0.3,win_minmax ); 
end

% Farneback method
of_level = 1; 
% of_nb = 5; 
% of_fs = 10; 
of_nb = 10;
of_fs = 20; 
opticFlow = opticalFlowFarneback( 'NumPyramidLevels',of_level, ... 
    'NeighborhoodSize',of_nb,'FilterSize',of_fs );
    
Vx = nan( size(I_bfld) ); 
Vy = nan( size(I_bfld) ); 
for t = 1:1:N
    flow_temp = estimateFlow( opticFlow, I_of(:,:,t) ); 
    Vx(:,:,t) = imgaussfilt( flow_temp.Vx, 1 ); 
    Vy(:,:,t) = imgaussfilt( flow_temp.Vy, 1 ); 
end

% Smooth the velocity field
% Vx = imgaussfilt3( Vx,[3 3 1] ); 
% Vy = imgaussfilt3( Vy,[3 3 1] ); 


% Coarse grain the velocity field
Vx_box = nan( length(Y),length(X),N ); 
Vy_box = nan( length(Y),length(X),N ); 
box_size = X(2)-X(1); 
for t = 1:1:N
    if rem( t,10 ) == 0
        disp( t )
    end
    Vx_box(:,:,t) = img_box( Vx(:,:,t),X,Y,box_size ); 
    Vy_box(:,:,t) = img_box( Vy(:,:,t),X,Y,box_size ); 
end

% Vx_box = imgaussfilt3( Vx_box,[3 3 1] ); 
% Vy_box = imgaussfilt3( Vy_box,[3 3 1] ); 

end


%% Local functions
% Coarse graining using mean
function I_box = img_box( I_in,x,y,box_size )
    m = length(y); 
    n = length(x); 
    I_box = zeros(m,n); 

    bs = round(box_size/2); 

    for i = 1:1:m
        for j = 1:1:n
            % I_box(i,j) = mean( I_in(y(i)-bs:y(i)+bs,x(j)-bs:x(j)+bs),'all' ); 
            I_box(i,j) = trimmean( I_in(y(i)-bs:y(i)+bs,x(j)-bs:x(j)+bs),5,'all' ); 
            % [N,xe] = histcounts( I_in(y(i)-bs:y(i)+bs,x(j)-bs:x(j)+bs),20 ); 
            % xc = xe(1:end-1)+(xe(2)-xe(1))/2; 
            % I_box(i,j) = xc( round(mean(find(N == max(N)))) ); 
        end
    end
end

% Image processing method by Katie
function I_out = method_Katie( l1,winGauss )
    % l1 = l1./imgaussfilt(l1,64);
    l1 = l1./imgaussfilt( l1,winGauss );
    l1 = imsharpen( l1,'Amount',3,'Radius',3 );
    l1 = imgaussfilt( l1,1 );
    I_out = normalise( l1 ); 
end


