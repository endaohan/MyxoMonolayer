% Revised based on Katie's code
% fpath - the folder for BFld
% l - the image (single frame)
% t - index of the image
% Do not save data files to the hard drive

% function kc_dfield(fpath,t,n)
function dirf = kc_dfield_noSave( l )

    n = 10; 

    % Calculate and save the director field for the experiment in fpath at
    % time t with smoothing kernel sized n.
    % addpath(genpath('../Misc. functions'));
    % [~,~,~] = mkdir(fpath,'analysis');
    % [~,~,~] = mkdir(fpath,'analysis/dfield/'); 
    
    % l = laserdata(fpath,t); 
    % Sharpen, which is vital
    ls = imsharpen(l,'Radius',3,'Amount',3);
    
    % Create gradient kernels for x and y.
    dl = 1/2*[-1 0 1];
    
    dx = zeros(numel(dl),numel(dl));
    dx(end/2+1/2,:) = dl;
    dy = zeros(numel(dl),numel(dl));
    dy(:,end/2+1/2) = dl;

    % Calculate the x and y gradients of the laser image.
    ls = padarray(ls,[1,1],'replicate','both');
    
    lx = conv2(ls,dx,'same');
    ly = conv2(ls,dy,'same');
    
    lx = lx(2:end-1,2:end-1);
    ly = ly(2:end-1,2:end-1);

    % Calculate the hetian for each pixel of the laser image from
    % gradients.
    
    J11 = lx.^2;
    J12 = lx.*ly;
    J22 = ly.^2;

    % Smooth it.
    
    J11 = (imgaussfilt(J11,n));
    J12 = (imgaussfilt(J12,n));
    J22 = (imgaussfilt(J22,n));

    % Calculate the eigenvector for the lower eigenvalue.
    
    VX = -J12;
    VY = (J11-J22)/2+sqrt((J22-J11).^2/4+J12.^2);

    % Save eigenvector direction as the director field.
    
    dirf = atan(VY./VX);
    
    %{
    name = sprintf('%06d.bin',t-1);
    fID = fopen([fpath '/analysis/dfield/' name],'w');
    fwrite(fID,dirf,'float');
    fclose(fID);
    %}
end

