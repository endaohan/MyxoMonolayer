function [X, Y, V] = PIV_peakDetection(I)
%
% PEAK_DETECTION  - Find the peak(s) (= local maxima) of a matrix, usually
% the outcome of a correlation. For multiple peaks, it finds them in order
% of their heigt, starting with the highest.
%
% usage:
%
% [V, X, Y] = peak_detection(I)
%
% I = input matrix to be searched for peaks
%
% output:
%
% X = x-position (column) of the peak
% Y = y-position (row) of the peak
% V = height of the peak

% the size of the input matrix
rows = size(I,1);
columns = size(I,2);

% we pre-allocate the arrays with an upper limit for the number of peaks,
% which is the number of data points in the input matrix I
X = zeros(rows*columns,1);
Y = zeros(rows*columns,1);
V = zeros(rows*columns,1);

% n is our counter for the number of peaks we find
n = 0;
for i=2:1:rows-1
    for j=2:1:columns-1
        
        % The following method is about 100x slower:
        % if I(i,j)==max( [I(j-1,i-1:i+1) I(j,i-1:i+1) I(j+1,i-1:i+1)] )
        
        % This is the 'fast' method. It looks a bit weird, but this 
        % if-statement is really the fastest method I found:
        if ...
                ( I(i,j) > I(i-1,j-1) ) && ...
                ( I(i,j) > I(i-1,j) ) && ...
                ( I(i,j) > I(i-1,j+1) ) && ...
                ( I(i,j) > I(i,j-1) ) && ...
                ( I(i,j) > I(i,j+1) ) && ...
                ( I(i,j) > I(i+1,j-1) ) && ...
                ( I(i,j) > I(i+1,j) ) && ...
                ( I(i,j) > I(i+1,j+1) )
            % we found a peak -> save position and value
            n = n + 1;
            X(n) = j;
            Y(n) = i;
            V(n) = I(i,j);
        end
    end
end


if n 
    % we found n peaks, remove rest of data from arrays
    X = X(1:n);
    Y = Y(1:n);
    V = V(1:n);
    % sort the peaks according to their height
    [V, IX] = sort(V,'descend');
    X = X(IX);
    Y = Y(IX);
else
    % we found no peak, return empty matrices
    X = [];
    Y = [];
    V = [];
end


