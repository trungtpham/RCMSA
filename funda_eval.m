clear
close all

%-------------------------------------------------------------------------%
model_type = 'fundamental8';
%-------------------------------------------------------------------------%

data_path = 'data/AdelaideRMF/F';

data_files = dir(data_path);
data_files(1:2) = [];
num_repetitions = 100;
seg_errors = zeros(length(data_files), num_repetitions);


for f=1:length(data_files)
    fprintf('Processing data %s \n', data_files(f).name);
    file_name = strcat(data_files(f).folder, '/', data_files(f).name);
    pack = load(file_name);
    
    %--- Prepare data --------------------------------------------------------%
    xy = pack.data;
    I1 = pack.img1;
    I2 = pack.img2;
    GT = pack.label;
    
    
    %---Normalize data--------------------------------------------------------%
    [dat_img_1, T1] = normalise2dpts(xy(1:3,:));
    [dat_img_2, T2] = normalise2dpts(xy(4:6,:));
    data = [dat_img_1 ; dat_img_2];
    %-------------------------------------------------------------------------%
    
    
    %----------Set parameters-------------------------------------------------%
    param.sig = 0.001;            % Standard deviation of noise
    param.smoothness = 0.25;      % Smoothness cost
    param.min_inliers = 10;       % Minimum number of inlier per structure
    param.rcm_sampling = 1;       % Used RCM sampling method
    param.sa    = 0.99;           % Simulated Annealing Schedule
    param.max_iteration = 5000;   % Max number of iterations
    param.min_iteration = 1000;   % Min number of iterations
    param.K     = 20;             % Patch size to update the weight
    %-------------------------------------------------------------------------%
    
      %---Robust model fitting--------------------------------------------------%
    
    for i=1:num_repetitions
        [estimated_pars, segmentation, energy] = rcmsa_model_fitting(data, xy, model_type, param);
        % Evaluation
        seg_errors(f, i) = segmentation_error(segmentation, GT);
    end
    fprintf('Median error = %f \n', median(seg_errors(f,:)));
    %-------------------------------------------------------------------------%
    
    %--Display segmentation result--------------------------------------------%
    display = 0;
    if display == 1
        figure(f);
        imshow(I1);hold on
        gscatter(xy(1,:), xy(2,:), segmentation, [], [], 20);
        title('The first label in red is outlier label');
    end
    drawnow;
    %-------------------------------------------------------------------------%
    
end

mean_err = mean(seg_errors, 2);
median_err = median(seg_errors, 2);
fprintf('\n');
fprintf('%s \t %s \t %s \n', pad('Image', 50), 'Mean Error', 'Median Error');
for f=1:length(data_files)
    fprintf('%s \t %f \t %f \n', pad(data_files(f).name, 50), mean_err(f), median_err(f));
end
