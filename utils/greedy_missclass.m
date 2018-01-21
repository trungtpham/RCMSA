function [total_miss] = greedy_missclass(Segmentation,npoints,ngroups)

% [miss,index] = missclass(Segmentation,npoints,ngroups)
%
% Computes the number of missclassified points in the vector Segmentation. 
%
% Segmentation: 1 by sum(npoints) or sum(ngroups) by 1 vector containing 
% the label for each group, ranging from 1 to n

% npoints: 1 by ngroups or ngroups by 1 vector containing the number of 
% points in each group.

% ngroups: number of groups
label_list = 1:ngroups;
npoints = [0 npoints];

for i=1:ngroups
    for k=1:length(label_list)
        label(k) = sum(Segmentation(sum(npoints(1:i))+1:sum(npoints(1:i+1)))==label_list(k));
    end
    [label_count ind] = max(label);
    label_list(ind)=[];
    miss(i) = npoints(i) - label_count;
    clear label
end

total_miss = sum(miss);

