%*************************************************************************%
% [GF, kx, ky] = FUN_GenerateGreenFun(Parameter, folder, lambda)
% 
% Given an MxN matrix where both M and N are even numbers, create the
% Green's function in Fourier space. 
% This new version allows changing the range of view of the Green's
% function. 
% Written by Endao Han, V8, 2021/1/16
% This version works for TFM experimental data analysis. 
%*************************************************************************%
function [GF, kx, ky] = FUN_GenerateGreenFun(Parameter, folder, lambda)

%% For debugging
%{
clear

% Route out
folder = 'E:\OneDrive - Princeton University\Projects\2018_FruitingBody\DATA\Confocal\20191030'; 

% Length, unit: m/pixel
mag = 1.5*60;       % Magnification of the microscope
Parameter.lscl = Param_lengthScale(mag)/1E6; 
% grid size
Parameter.d = 1; 
Parameter.box_expt = Parameter.lscl*Parameter.d; 

% Parameters that match our experiments with PAA
% Young's modulus
Parameter.G = 230; 
% Poisson's ratio
Parameter.nu = 0.5; 
Parameter.E = 2*( 1+Parameter.nu )*Parameter.G; 

% Size of the system
Ng = 200; 
Parameter.N_grid = [2*Ng+1 2*Ng+1]; 

% Noise term
lambda = 0; %1/1E10; 
%}



%% Basic parameters
N_grid = Parameter.N_box; 
if ~exist('lambda','var')
    lambda = 0; 
end
Lscl = Parameter.box_expt; 
nu = Parameter.nu; 
E = Parameter.E; 

fileName = 'GreensFunction'; 
routeOut = [folder '\' fileName]; 
mkdir(routeOut)

% Check if the dimensions are odd numbers
if ~( rem(N_grid(1),2) && rem(N_grid(2),2) )
    disp('ERROR: FUN_GenerateGreenFun: matrix dimension is not ODD! ')
    return
end

% Save results or not: 1 - save; 0 - do not save
flagSave = 1; 



%% Calculate Green's functions
%{
% Point force in x
F_x = zeros(N_grid(1),N_grid(2),2); 
F_x(ceil(N_grid(1)/2),ceil(N_grid(2)/2),:) = [1 0];        % unite: N

OUT_x = FUN_ForceToDisp_ODD(Parameter,F_x,[folder '\Disp_Fx1N'],0); 

G_11 = fftshift(fft2(OUT_x.U_m)); 
G_21 = fftshift(fft2(OUT_x.V_m)); 

% Point force in y
F_y = zeros(N_grid(1),N_grid(2),2); 
F_y(ceil(N_grid(1)/2),ceil(N_grid(2)/2),:) = [0 1];        % unite: N

OUT_y = FUN_ForceToDisp_ODD(Parameter,F_y,[folder '\Disp_Fy1N'],0); 

G_12 = fftshift(fft2(OUT_y.U_m)); 
G_22 = fftshift(fft2(OUT_y.V_m)); 
%}

% Check the z position of the slice
if ~exist('Parameter.Z','var')
    Z = 0; 
else
    Z = Parameter.Z; 
end

% Calculate k vectors
Ng = (N_grid-1)/2; 
ky = linspace( -Ng(1),Ng(1),N_grid(1) )*2*pi/Lscl/N_grid(1); 
kx = linspace( -Ng(2),Ng(2),N_grid(2) )*2*pi/Lscl/N_grid(2); 
[Kx, Ky] = meshgrid(kx,ky); 

% Green's function in Fourier space
G_pre = 2*(1+nu)*exp(Z*sqrt(Kx.^2+Ky.^2))/E./(Kx.^2+Ky.^2).^1.5; 
G_11 = G_pre.*( Kx.^2 + Ky.^2 - Kx.*Kx.*(nu-Z*sqrt(Kx.^2+Ky.^2)/2) ); 
G_12 = G_pre.*( - Kx.*Ky.*(nu-Z*sqrt(Kx.^2+Ky.^2)/2) ); 
G_21 = G_12; 
% G_21 = G_pre.*( - Kx.*Ky.*(nu-Z*sqrt(Kx.^2+Ky.^2)/2) ); 
G_22 = G_pre.*( Kx.^2 + Ky.^2 - Ky.*Ky.*(nu-Z*sqrt(Kx.^2+Ky.^2)/2) ); 

G = cell(2,2); 
G{1,1} = G_11; 
G{1,2} = G_12; 
G{2,1} = G_21; 
G{2,2} = G_22; 


% M matrix
LbdSqr = ones( N_grid(1),N_grid(2) )*lambda^2; 

Mat.a = G_11.^2+G_21.^2+LbdSqr; 
Mat.b = G_11.*G_12+G_21.*G_22; 
Mat.c = Mat.b; 
Mat.d = G_12.^2+G_22.^2+LbdSqr; 

Denom = ( Mat.a.*Mat.d - Mat.b.*Mat.c ).^(-1); 

M_11 = Denom.*(Mat.d.*G_11-Mat.b.*G_12); 
M_12 = Denom.*(Mat.d.*G_21-Mat.b.*G_22); 
M_21 = Denom.*(-Mat.c.*G_11+Mat.a.*G_12); 
M_22 = Denom.*(-Mat.c.*G_21+Mat.a.*G_22); 

% Set the central values of the M matrix to be zero
M_11( Ng(1)+1, Ng(2)+1 ) = 0; 
M_12( Ng(1)+1, Ng(2)+1 ) = 0; 
M_21( Ng(1)+1, Ng(2)+1 ) = 0; 
M_22( Ng(1)+1, Ng(2)+1 ) = 0; 

M = cell(2,2); 
M{1,1} = M_11; 
M{1,2} = M_12; 
M{2,1} = M_21; 
M{2,2} = M_22; 


% Output
GF.G = G; 
GF.M = M; 



%% Present the results

% Plot Green's function in Fourier space
hd_GF = figure; 
subplot(2,2,1)
imagesc( kx,ky, G{1,1} )
ax = gca;
ax.YDir = 'normal'; 
title( '$\tilde{G}_{11}$','Interpreter','Latex' )
xlabel('k_x (m^{-1})')
ylabel('k_y (m^{-1})')
axis equal
axis([min(kx),max(kx),min(ky),max(ky)])
colorbar

subplot(2,2,2)
imagesc( kx,ky,G{1,2} )
ax = gca;
ax.YDir = 'normal'; 
title( '$\tilde{G}_{12}$','Interpreter','Latex' )
xlabel('k_x (m^{-1})')
ylabel('k_y (m^{-1})')
axis equal
axis([min(kx),max(kx),min(ky),max(ky)])
colorbar

subplot(2,2,3) 
imagesc( kx,ky, G{2,1} )
ax = gca;
ax.YDir = 'normal'; 
title( '$\tilde{G}_{21}$','Interpreter','Latex' )
xlabel('k_x (m^{-1})')
ylabel('k_y (m^{-1})')
axis equal
axis([min(kx),max(kx),min(ky),max(ky)])
colorbar

subplot(2,2,4)
imagesc( kx,ky, G{2,2} )
ax = gca;
ax.YDir = 'normal'; 
title( '$\tilde{G}_{22}$','Interpreter','Latex' )
xlabel('k_x (m^{-1})')
ylabel('k_y (m^{-1})')
axis equal
axis([min(kx),max(kx),min(ky),max(ky)])
colorbar


% Plot the inversed Green's function M
hd_M = figure; 
subplot(2,2,1)
imagesc( kx,ky, M{1,1} )
ax = gca;
ax.YDir = 'normal'; 
title( '$\tilde{M}_{11}$','Interpreter','Latex' )
xlabel('k_x (m^{-1})')
ylabel('k_y (m^{-1})')
axis equal
axis([min(kx),max(kx),min(ky),max(ky)])
colorbar

subplot(2,2,2)
imagesc( kx,ky,M{1,2} )
ax = gca;
ax.YDir = 'normal'; 
title( '$\tilde{M}_{12}$','Interpreter','Latex' )
xlabel('k_x (m^{-1})')
ylabel('k_y (m^{-1})')
axis equal
axis([min(kx),max(kx),min(ky),max(ky)])
colorbar

subplot(2,2,3) 
imagesc( kx,ky, M{2,1} )
ax = gca;
ax.YDir = 'normal'; 
title( '$\tilde{M}_{21}$','Interpreter','Latex' )
xlabel('k_x (m^{-1})')
ylabel('k_y (m^{-1})')
axis equal
axis([min(kx),max(kx),min(ky),max(ky)])
colorbar

subplot(2,2,4)
imagesc( kx,ky, M{2,2} )
ax = gca;
ax.YDir = 'normal'; 
title( '$\tilde{M}_{22}$','Interpreter','Latex' )
xlabel('k_x (m^{-1})')
ylabel('k_y (m^{-1})')
axis equal
axis([min(kx),max(kx),min(ky),max(ky)])
colorbar



%% Save results
if flagSave == 1
    fileName = ['_lambda' num2str(lambda) '_boxSize' num2str(Parameter.d)]; 
    saveas(hd_GF, [routeOut '_G_' fileName '.fig'])
    saveas(hd_GF, [routeOut '_G_' fileName '.png'])

    saveas(hd_M, [routeOut '_M_' fileName '.fig'])
    saveas(hd_M, [routeOut '_M_' fileName '.png'])

    save([routeOut '_' fileName '.mat'],'GF','kx','ky')

    close(hd_GF)
    close(hd_M)
end






