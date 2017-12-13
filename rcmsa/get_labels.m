function [label, eng, pcost] = get_labels(graph, datacost, smoothcost)
% Input
% graph [NxN]     : represents a sparse adjacency matrix
% datacost [NxL]  : contains cost of N points w.r.t L labels
% datacost_scale  : scale of data cost
% smoothcost      : smoothness cost of the MRF model 

% Output
% label           : the output labels
% eng             : the energy of the system with the esimated labels
% pcost           : cost of each point w.r.t its assigned label.


% energy scale required by GCO
energy_scale = 1e3;

% Make squared residual
datacost = datacost.^2;

smoothcost = smoothcost.*energy_scale;
datacost = min(1e7, datacost.*energy_scale);

% Label initalisation
[n, num_hyp ] = size(datacost);
[~, ilabel] = min(datacost,[],2);

% Create Graph cut object
h = GCO_Create(n,num_hyp);

% Uniform cost is used for all label pairs
S = ~eye(num_hyp);
S(1,:) = 0.1;
S(:,1) = 0.1;
GCO_SetSmoothCost(h,int32(smoothcost.*S));

% Set neighbors
GCO_SetNeighbors(h,graph);

% Set data cost
GCO_SetDataCost(h,int32(datacost'));
%GCO_SetLabelCost(h,int32(0)); 
% Set labels
GCO_SetLabeling(h,ilabel)


% run alpha-expansion
GCO_Expansion(h);                

% Computer energy
[E, D, S, L] = GCO_ComputeEnergy(h);

% Output engery
eng = double(E)/energy_scale;

% get labelling results
label = GCO_GetLabeling(h);

% Get cost for each datum w.r.t to their label
Tdcost = datacost'; % Transpose dcost matrix
ind = dummyvar(double(label))';
pcost = Tdcost(logical(ind)); % cost of data points w.r.t their labels
pcost = pcost/energy_scale;
% Delete Graph Cut object
GCO_Delete(h);

end
