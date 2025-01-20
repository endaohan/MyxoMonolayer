function J = PIVBlueWhiteRed(m)

if nargin < 1
   m = size(get(gcf,'colormap'),1);
end

r = [linspace(0,1,ceil(m/2))'; ones(floor(m/2),1)];
b = [ones(floor(m/2),1); linspace(1,0,ceil(m/2))'];
g = [r(1:ceil(m/2)); b(ceil(m/2)+1:end)];

J = [r g b];
