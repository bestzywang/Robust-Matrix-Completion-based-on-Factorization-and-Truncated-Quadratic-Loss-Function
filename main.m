clear
clc
close all
peaksnr_1 = [];
peaksnr_2 = [];

ssim_1 = [];
ssim_2 = [];

t_1 = [];
t_2 = [];



for ii=1:100

image = imread('18.jpg');
[width,height,z]=size(image);
if(z>1)
    image=rgb2gray(image);
end
% unit8 to double
image = mat2gray(image);
[m,n] = size(image);
M = image;

% set the estimated rank value
rak = 10;

% Set the percentage of observed entries
per = 0.8;
array_Omega = binornd( 1, per, [ m, n ] );

% add mixture noise
M_noise1 = M + 0.02*randn(m,n);
M_noise = imnoise(M_noise1,"salt & pepper",0.1);
M_noise = M_noise.*array_Omega;
   

maxiter = 50;

% ip = 4 for real images, ip=3 for synthetic data.
ip = 4;

tic
[X_A,~,~,RMSE_0 ]= HOAT(M_noise,array_Omega,rak,maxiter,ip);
toc
t_1 = [t_1 toc];
peaksnr_1 = [peaksnr_1 psnr(M,X_A)];
ssim_1 = [ssim_1 ssim(M,X_A)];


%% HOMT

tic
[X_M,~,~,RMSE_1] = HOMT(M_noise,array_Omega,rak,maxiter,ip);
toc
t_2 = [t_2 toc];
peaksnr_2 = [peaksnr_2 psnr(M,X_M)];
ssim_2 = [ssim_2 ssim(M,X_M)];

end

