%******************************************************************%
%         2D JAMMING FRONT: How does the front get flat
%          James Frank Institute, University of Chicago
%                 Written by: Endao Han, 20130501
%******************************************************************%

function bin_ave_data = bin_ave(x,y,n_bin)

x_min = min(x);
x_max = max(x);
y_min = min(y);
y_max = max(y);

x_leng = length(x);
y_leng = length(y);
if x_leng == y_leng
    leng = x_leng;
else
    sprintf('invalid data')
end    

%************************************************************************%
% Pick out the particles on the edge two bins and in the middle two bins,
% and put them in two separate arrays.
%************************************************************************%

bin = (max(x)-min(x))/n_bin;
bin_ave_data = zeros(n_bin,5);

for i = 1:1:n_bin
    ind = find((x>=(i-1)*bin)&(x<i*bin));
    y_temp = y(ind);
    x_bin = (i-0.5)*bin;
    y_bin = mean(y_temp);
    std_bin = std(y_temp);
    num_bin = length(ind);
    
    bin_ave_data(i,1) = i;
    bin_ave_data(i,2) = x_bin;
    bin_ave_data(i,3) = y_bin;
    bin_ave_data(i,4) = std_bin;
    bin_ave_data(i,5) = num_bin;
end
