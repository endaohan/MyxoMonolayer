function S = kc_nemorderfield_noSave( dfield )
% Calculates local order field of image in fpath, at time t.
%    Saves into analysis/order folder as an array of floats.
%% Calc order
    %r = 0.5;
    r = 0.5;
    % XYcal = getXYcal(fpath);
    XYcal = 0.1; 
    %{
    name = sprintf('%06d.bin',t-1);
    [~,~,~] = mkdir([fpath '\analysis\order']);
    try
        dfield = kc_loaddata(fpath,t,'dfield','float',dImg);
    catch
        msgbox(['Run dfield calculations for ' fpath ' at time t = ' t '...']);
        return
    end
    %}
    rr = round(r/XYcal);
    se = strel('disk',rr);
    mf = se.Neighborhood;
    norm = conv2(ones(size(dfield)),mf,'same');

    dfield2 = 2*dfield;
    dx2 = cos(dfield2);
    dy2 = sin(dfield2);

    mdx2 = conv2(dx2,mf,'same')./norm;
    mdy2 = conv2(dy2,mf,'same')./norm;

    mdf = atan2(mdy2,mdx2)/2;

    %Prefered direction (n).
    nx = cos(mdf);
    ny = sin(mdf);

    %Director field mean over areas.
    dx = cos(dfield);
    dy = sin(dfield);

    mdxdx = conv2(dx.*dx,mf,'same')./norm;
    mdxdy = conv2(dx.*dy,mf,'same')./norm;
    mdydy = conv2(dy.*dy,mf,'same')./norm;

    S = (2*(nx.*nx.*mdxdx+2*nx.*ny.*mdxdy+ny.*ny.*mdydy)-1);

    %% Save file
    % fID = fopen([fpath '\analysis\order\' name],'w');
    % fwrite(fID,S,'float');
    % fclose(fID);
    
end