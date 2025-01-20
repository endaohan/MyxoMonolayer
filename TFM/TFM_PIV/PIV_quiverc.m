 function hh = PIV_quiverc(x, y, u, v, cmap, scale)
% Modified version of Quiver to plots velocity vectors as arrows 
% with components (u,v) at the points (x,y) using the current colormap 


%PIV_quiverc Quiver color plot.
%   PIV_quiverc(x, y, u, v, cmap, scale) plots velocity vectors as arrows 
%   with components (u,v) at the points (x,y).  
%   The matrices x,y,u,v must all be the same size.
%   cmap determines the color that will be given to the vectors.
%   Vectors are stretched by an amount scale

%-------------------------------------------------------------
% Arrow head parameters
alpha = 0.33; % Size of arrow head relative to the length of the vector
beta = 0.23;  % Width of the base of the arrow head relative to the length

%----------------------------------------------
% Define colormap
vr=sqrt(u.^2+v.^2);
% minimum velocity
minV = min(vr(:));
% maximum velocity
maxV = max(vr(:));
vr = vr - minV;
vr = vr / (maxV-minV);
vrn = uint8(vr * 32);
% vrn=round(vr/max(vr(:))*64);
CC=cmap;
ax = newplot;
next = lower(get(ax,'NextPlot'));
hold_state = ishold;

%----------------------------------------------
% Scale the velocities
u = u * scale;
v = v * scale;

%----------------------------------------------
% Make velocity vectors and plot them

x = x(:).';y = y(:).';
u = u(:).';v = v(:).';
vrn=vrn(:).';
uu = [x; x+u; NaN(size(u))];
vv = [y; y+v; NaN(size(u))];
vrn1= [vrn; NaN(size(u)); NaN(size(u))];

uui=uu(:);  vvi=vv(:);  vrn1=vrn1(:); imax=size(uui);
hold on

for i=1:3:imax-1
    ii=uint8(round(vrn1(i)));
    if ii==0; 
        ii=1; 
    end        
    c1= CC(ii,1);
    c2= CC(ii,2);
    c3= CC(ii,3);
    plot(uui(i:i+1),vvi(i:i+1),'color',[c1 c2 c3]);
end

%----------------------------------------------
% Make arrow heads and plot them
hu = [x+u-alpha*(u+beta*(v+eps));x+u; ...
    x+u-alpha*(u-beta*(v+eps));NaN(size(u))];
hv = [y+v-alpha*(v-beta*(u+eps));y+v; ...
    y+v-alpha*(v+beta*(u+eps));NaN(size(v))];
vrn2= [vrn;vrn;vrn;vrn];

uui=hu(:);  vvi=hv(:);  vrn2=vrn2(:); imax=size(uui);

for i=1:imax-1
    ii=uint8(round(vrn2(i)));
    if ii==0
        ii=1;
    end
    c1= CC(ii,1);
    c2= CC(ii,2);
    c3= CC(ii,3);
    plot(uui(i:i+1),vvi(i:i+1),'color',[c1 c2 c3]);
end

%----------------------------------------------

if ~hold_state
    hold off
    view(2)
    set(ax,'NextPlot',next)
end

if nargout>0
    hh = [h1;h2;h3];
end