function err = segmentation_error(pred, groundtruth)

if min(pred) == 0
    pred = pred + 1;
end

if min(groundtruth) == 0
    groundtruth = groundtruth + 1;
end

ulabels = unique(groundtruth);
gtrue_num_models = length(ulabels);
idx = [];
for l=1:gtrue_num_models
    idx = [idx, find(groundtruth == ulabels(l))];
end

pred = pred(idx);
groundtruth = groundtruth(idx);

num_per_model = histc(groundtruth, unique(groundtruth));

if gtrue_num_models<10
    mc = missclass(pred, num_per_model, gtrue_num_models);
    mc = mc/length(pred)*100;
else
    mc = greedy_missclass(pred, num_per_model, gtrue_num_models);
    mc = mc/length(pred)*100;
end
err = mc;
end