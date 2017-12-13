function [V_R]  = rcm_sampling(data, psize, degenfn, f, e, q, cost)
% Follows Swendsen-Wang cut of Barbu and Zhu.
%-------
% Input:
%-------
% f (Nx1)       : Current labelling.
% e (|E|x2)     : Edges.
% q (|E|x1)     : Edge probabilities.
% psize         : minimal number of points required to construct a model
% degenfn       : degerate function
%--------
% Output:
%--------
% V_R    : Vertices in R.

% Constants.
N = length(f);

% Sample the connected component.
cost = sqrt(cost);

% Edges that remain on due to same labels.
eon_det = f(e(:,1))==f(e(:,2));

% Edges that are turned on stochastically.
eon_sto = rand(length(q),1)<=q;

% Either the edge is already off due to different labels, or
% the edge is turned off stochastically.
eon = logical(eon_det.*eon_sto);

% Get current set of connected components.
Eon = sparse(e(eon,1),e(eon,2),ones(sum(eon),1),N,N);
Eon = Eon + Eon';   % Make symmetric.
[S,C] = graphconncomp(Eon);

% Compute sampling weights
w = ones(1,S);
for s=1:S
    w(s) = sum(cost(C==s));
end
counts = histc(C, unique(C));
w = w./counts;
w(counts<psize) = 0;

pick_valid = 0;
max_num_trials = 1000;
trials = 0;
while pick_valid == 0
    trials = trials + 1;
    if trials > max_num_trials
        break;
    end
    
    % Pick a connected component R probabilistically.
    if sum(w) > 0
        R = randsample(S,1,true,w);
        V_R = find(C==R);
    else
        V_R = randsample(N, psize);
    end
       
    % While sampling based on the current fitting error is more effective, the
    % cost of computing sampling weights is a bit expensive. 
    % Using the below line if wanna sample randomly
    % R = randsample(S,1);
    
    
    
    if length(V_R) < psize
        continue;
    end
    
    % Check validation of the sampled subset
    isdegen= feval(degenfn,data(:,V_R));
    if isdegen==0
        pick_valid = 1;
    end
end


