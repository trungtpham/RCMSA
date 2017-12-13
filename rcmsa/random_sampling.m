function [V_R]  = random_sampling(data, psize, degenfn)
% psize         : minimal number of points required to construct a model
% degenfn       : degerate function
%--------
% Output:
%--------
% V_R    : Vertices in R.

N = length(data);
pick_valid = 0;
while pick_valid == 0
    
    V_R = randsample(N, psize);
    
    if length(V_R) < psize
        continue;
    end    
    % Check validation of the sampled subset
    isdegen= feval(degenfn,data(:,V_R));
    if isdegen==0
        pick_valid = 1;
    end
end


