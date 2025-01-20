% ************************************************************************
% disp2traction( u_m,u_m,GF )
% Input: u_m - displacement field in x
% v_m - displacement field in y
% GF - Green's function (structure)
% Output: [f_x, f_y] in the unit of Pa
% ************************************************************************

function [f_x,f_y] = disp2traction( u_m, v_m, GF )

% Check if the matrices have the same size
Du = size(u_m); 
Dv = size(v_m); 
if Du~=Dv
    disp("ERROR: disp2traction, input matrices not the same size! ")
    return
end
% Check if the dimensions are both odd
if sum(rem(Du,2)) ~= 2
    disp("ERROR: disp2traction, input matrices is not odd! ")
    return
end

% Transformation
uFFT = fftshift( fft2(u_m) );
vFFT = fftshift( fft2(v_m) ); 

fFFT_x = GF.M{1,1}.*uFFT+GF.M{1,2}.*vFFT; 
fFFT_y = GF.M{2,1}.*uFFT+GF.M{2,2}.*vFFT; 

f_x = ifft2( ifftshift(fFFT_x),'symmetric' ); 
f_y = ifft2( ifftshift(fFFT_y),'symmetric' ); 

end
    
    