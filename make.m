%--Add path---------------------------------------------------------------%
addpath('model_specific');
addpath('gco-v3.0');
addpath('rcmsa');
addpath('gco-v3.0/matlab');
%-------------------------------------------------------------------------%

%-------------Compile required mex file-----------------------------------%
if (exist('computeIntersection')~=3)
    mex rcmsa/computeIntersection.c
end
%-------------------------------------------------------------------------%
fprintf('Done!\n')
