function [ V_R]  = rcm_sampling(data, psize, degenfn, f, e, q, cost)
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
V_R = [];
valid = 0;
while valid == 0
    
    
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
    
    % Pick a connected component R probabilistically.
    w = ones(1,S);
    for s=1:S
        if sum(C==s) < 8 % Remove small clusters having less than 8 points
            w(s) = 1e-100;
        else 
            w(s) = mean(cost(C==s));
        end
    end
    
    R = randsample(S,1,true,w);
    % Using the below line if wanna sample randomly
    % R = randsample(S,1);
    
    V_R = find(C==R);
    
    % Check validation of the sampled subset
    isdegen= feval(degenfn,data(:,V_R));
    if length(V_R) > psize && isdegen==0
        valid = 1;
    end
end


