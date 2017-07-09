function [ fitfn resfn degenfn psize numpar ] = getModelPara(model_type)

%---------------------------
% Model specific parameters.
%---------------------------

switch model_type
    
    case 'line'
        fitfn = @line_fit;
        resfn = @line_res;
        degenfn = @line_degen;
        psize = 2;
        numpar = 3;
    case 'circle'
        fitfn = @circle_fit;
        resfn = @circle_res;
        degenfn = @circle_degen;
        psize = 3;
        numpar = 3;
    case 'homography'
        fitfn = @homography_fit;
        resfn = @homography_res;
        degenfn = @homography_degen;
        psize = 4;
        numpar = 9;
    case 'fundamental'
        fitfn = @fundamental_fit;
        resfn = @fundamental_res;
        degenfn = @fundamental_degen;
        psize = 7;
        numpar = 9;
    case 'fundamental8'
        fitfn = @fundamental_fit8;
        resfn = @fundamental_res;
        degenfn = @fundamental_degen;
        psize = 8;
        numpar = 9;
    case 'motion'
        fitfn = @motion_fit;
        resfn = @motion_res;
        degenfn = @motion_degen;
        psize = 4;
        numpar = 25;
    otherwise
        error('unknown model type!');
end

end