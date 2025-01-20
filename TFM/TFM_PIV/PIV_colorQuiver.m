function PIV_colorQuiver(PIVParams, Vectors)

% the number of different colors that we want to show
% 32 is really enough for color-coding vectors
numSteps = 32;

% get the data
x = Vectors.x;
y = Vectors.y;
u = Vectors.u;
v = Vectors.v;
u = u * PIVParams.VectorScale;
v = v * PIVParams.VectorScale;

% the value that is going to determine the colors
z = Vectors.magnitude;

% find vectors that we shouldn't plot:
validIndices = logical(Vectors.validationFlag~=1);
% do not plot invalid vectors:
x = x(validIndices);
y = y(validIndices);
u = u(validIndices);
v = v(validIndices);
z = z(validIndices);

% the range that the different colors will match:
% minZ = min(z(:));
% maxZ = max(z(:));
minZ = PIVParams.ScalarMinValue;
maxZ = PIVParams.ScalarMaxValue;
display('Color coded vectors.')
display(['Minimum value: ' num2str(min(min(z)))])
display(['Maximum value: ' num2str(max(max(z)))])
display(['Mean value:    ' num2str(mean(mean(z)))])

% determine the step size
zStep = (maxZ - minZ) / numSteps;

% select a colormap
switch PIVParams.VectorColormap
    case 'BlueWhiteRed'
        load PIVColors_BlueWhiteRed.mat
    case 'GreenYellowRed'
        load PIVColors_GreenYellowRed.mat
    otherwise
        figure
        cmap = eval(['colormap(' PIVParams.VectorColormap '(numSteps))']);
        close
end

% % interpolate the colormap to match numSteps:
% cmapNew(:,1) = interp1(cmap(:,1),(1:256)',xi');
% cmapNew(:,2) = interp1(cmap(:,2),(1:256)',xi');
% cmapNew(:,3) = interp1(cmap(:,3),(1:256)',xi');
% cmap = cmapNew;

% bin the vectors
for i=2:numSteps-1
    bin = z>=(minZ+(i-1)*zStep) & z<(minZ+i*zStep);
    quiver(x(bin), y(bin), u(bin), -v(bin), 0, 'color', cmap(i,:))
end
bin = z>(minZ+(numSteps-1)*zStep);% & z<=(minZ+numSteps*zStep);
quiver(x(bin), y(bin), u(bin), -v(bin), 0, 'color', cmap(i,:))