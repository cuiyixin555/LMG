clc;
clear;
close all;
addpath(genpath('whyte_code'));
addpath(genpath('cho_code'));
addpath(genpath('chen code'));
opts.prescale = 1; %%downsampling
opts.xk_iter = 5; %% the iterations
opts.gamma_correct = 1.0;
opts.k_thresh = 20;

filename = './test_data/im01_ker01.mat'; 
opts.kernel_size = 27;  saturation = 0;
lambda_lmg =4e-3; lambda_grad =4e-3;opts.gamma_correct = 1;
lambda_tv = 0.001; lambda_l0 = 5e-4; weight_ring = 1;
load(filename)
%%
%===================================

% y = imread(filename);
% y = y(3:end-2,3:end-2,:);
% y = imfilter(y,fspecial('gaussian',5,1),'same','replicate'); 
isselect = 0; %false or true
if isselect ==1
    figure, imshow(y);
    %tips = msgbox('Please choose the area for deblurring:');
    fprintf('Please choose the area for deblurring:\n');
    h = imrect;
    position = wait(h);
    close;
    B_patch = imcrop(y,position);
    y = (B_patch);
else
    y = y;
end
if size(y,3)==3
    yg = im2double(rgb2gray(y));
else
    yg = im2double(y);
end
tic;
[kernel, interim_latent] = blind_deconv(yg, lambda_lmg, lambda_grad, opts);
toc
%% Algorithm is done!
%% ============Non-blind deconvolution ((Just use text deblurring methods))====================%%
y = im2double(y);
%% Final Deblur: 
if ~saturation
    Latent = ringing_artifacts_removal(y, kernel, 2e-3, 1e-3, 0);
else
    %% 2. Whyte's deconvolution method (For saturated images)
    Latent = whyte_deconv(y, kernel);
end
figure; imshow(Latent)
imwrite(Latent,'chen.png');
k = kernel/max(kernel(:));
imwrite(k,'chengk.jpg');

