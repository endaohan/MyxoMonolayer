function Vectors = PIV_calcScalarFields(PIVParams, Vectors)

%
% PIV_calcScalarFields - calculate scalars from the velocity field. The
% scalars that will be calculated are given by the parameters file.
% (This function is called by the function PIV)
%


% read and copy vector data
x = Vectors.x;
y = Vectors.y;
u = Vectors.u;
v = Vectors.v;

% get dimensions of vector field
width = find(y==y(1),1,'last');
% width = i-1;
height = length(x) / width;

% reshape vectors into matrices:
% x = reshape(x,width,height)';
% y = reshape(y,width,height)';
u = reshape(u,width,height)';
v = reshape(v,width,height)';

% get horizontal and vertical spacing
dx = PIVParams.GridSize - PIVParams.Overlap;
dy = PIVParams.GridSize - PIVParams.Overlap;

% get correct scaling for the spacing
dx = dx * PIVParams.Scaling;
dy = dy * PIVParams.Scaling;

% calculate velocity gradients
%     [dudx, dudy] = gradient(u,dx,dy);
%     [dvdx, dvdy] = gradient(v,dx,dy);
% use custom functions to calculate gradients using least squares
% approach. This is a suitable method for PIV-data
% See PIV a practical guide, second edition, Springer
% differentiation filters (x- and y-direction):
fx = [-2 -1 0 1 2];
fy = [-2; -1; 0; 1; 2];
dudx = imfilter(u,fx,'replicate') / (10*dx);
dudy = imfilter(u,fy,'replicate') / (10*dy);
dvdx = imfilter(v,fx,'replicate') / (10*dx);
dvdy = -imfilter(v,fy,'replicate') / (10*dy);



% set values along the edges to zero because we only trust the central
% differences
% we do this with a mask:
mask = zeros(height,width);
mask(2:end-1,2:end-1) = 1;
dudx = dudx .* mask;
dudy = dudy .* mask;
dvdx = dvdx .* mask;
dvdy = dvdy .* mask;

% calculate velocity magnitude
magnitude = sqrt(u.^2 + v.^2);

% calculate vorticity
vorticity = dudy - dvdx;
shearStrain = dudy + dvdx;
normalStrain = dudx + dvdy;

% reshape matrices to vectors to be compatible to the other data
% number of vectors N:
N = width * height;
Vectors.magnitude = reshape(magnitude',1,N);
Vectors.signalToNoiseRatio = zeros(1,N);
Vectors.vorticity = reshape(vorticity',1,N);
Vectors.shearStrain = reshape(shearStrain',1,N);
Vectors.normalStrain = reshape(normalStrain',1,N);
Vectors.dudx = reshape(dudx',1,N);
Vectors.dudy = reshape(dudy',1,N);
Vectors.dvdx = reshape(dvdx',1,N);
Vectors.dvdy = reshape(dvdy',1,N);
