% Fundamental matrix fitting code adapted from Peter Kovesi's
% implementation of 7-point fundamental matrix estimation. This code makes
% use of Andrew Zisserman's 7 point fundamental matrix code.
% See:  http://www.robots.ox.ac.uk/~vgg/hzbook/code/

function P = fundamental_fit(X)

Fvgg = vgg_F_from_7pts_2img(X(1:3,:), X(4:6,:));

if isempty(Fvgg)
    P = [];
    return;
end

% Store the (potentially) 3 solutions in a cell array
[rows,cols,Nsolutions] = size(Fvgg);

P = cell(Nsolutions,1);
for n = 1:Nsolutions
    P{n} = Fvgg(:,:,n);
end

end


    