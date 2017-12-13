clear
close all
%clc

%-------------------------------------------------------------------------%
model_type = 'homography';
%-------------------------------------------------------------------------%

%--Load data--------------------------------------------------------------%
pack = load('./data/raglan.mat');
%pack = load('./data/unionhouse.mat');
%pack = load('./data/mc1.mat');
%pack = load('./data/mc3.mat');
%-------------------------------------------------------------------------%

%--- Prepare data --------------------------------------------------------%
xy = pack.data;
I1 = pack.img1;
I2 = pack.img2;

%---Normalize data--------------------------------------------------------%
[dat_img_1, T1] = normalise2dpts(xy(1:3,:));
[dat_img_2, T2] = normalise2dpts(xy(4:6,:));
data = [ dat_img_1 ; dat_img_2 ];
%-------------------------------------------------------------------------%

%----------Set parameters-------------------------------------------------%
param.sig = 0.03;             % Standard deviation of noise
param.min_inliers = 25;       % Minimum number of inlier per structure
param.rcm_sampling = 1;       % Used RCM sampling method
param.sa    = 0.9;            % Simulated Annealing Schedule
param.M     = 2000;           % Max number of iterations
param.K     = 100;             % Patch size to update the weight
%-------------------------------------------------------------------------%

%---Robust model fitting--------------------------------------------------%
[estimated_pars, segmentation] = rcmsa_model_fitting(data, xy, model_type, param);
%-------------------------------------------------------------------------%


%--Display ground truth segmentation ---------------------------------%
display = 1;
%if display == 1
%    figure(1);
%    imshow(I1);hold on
%    gscatter(xy(1,:), xy(2,:), GT, [], [], 20);
%end
%-------------------------------------------------------------------------%

%--Display segmentation result--------------------------------------------%
if display == 1
    figure(1);
    imshow(I1);hold on
    gscatter(xy(1,:), xy(2,:), segmentation, [], [], 20);
    title('The first label in red is outlier label');
end
drawnow;
%-------------------------------------------------------------------------%
