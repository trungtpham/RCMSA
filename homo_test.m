clear 
close all
clc

%-------------------------------------------------------------------------%
model_type = 'homography';
%-------------------------------------------------------------------------%

%--Load data--------------------------------------------------------------%
%pack = load('./data/raglan.mat');
%pack = load('./data/unionhouse.mat');
%pack = load('./data/mc1.mat');
pack = load('./data/mc3.mat');
%-------------------------------------------------------------------------%

%--- Prepare data --------------------------------------------------------%
xy = pack.data;
I1 = pack.img1;
I2 = pack.img2;

%---Normalize data--------------------------------------------------------%
[dat_img_1 T1] = normalise2dpts(xy(1:3,:));
[dat_img_2 T2] = normalise2dpts(xy(4:6,:));
data = [ dat_img_1 ; dat_img_2 ];
%-------------------------------------------------------------------------%


%----------Set parameters-------------------------------------------------%
param.sig = 0.005;    % Standard deviation of noise (inlier threshold)
param.bet = 10;       % Model complextiy penalty
param.sa  = 0.999;    % Simulated Annealing Schedule
param.M   = 1000;     % Max number of iterations
param.K   = 20;       % Patch size to update the weight
%-------------------------------------------------------------------------%

%---Robust model fitting--------------------------------------------------%
[epar elabel] = rcmsa_model_fitting(data, xy, model_type, param);         
%-------------------------------------------------------------------------%

%--Display segmentation result--------------------------------------------%
figure;
imshow(I1);hold on
gscatter(xy(1,:), xy(2,:), elabel, [], [], 20);
title('The first label in red is outlier label');
figure(2);
imshow(I2);hold on
gscatter(xy(4,:), xy(5,:), elabel, [], [], 20);
title('The first label in red is outlier label');
%-------------------------------------------------------------------------%
