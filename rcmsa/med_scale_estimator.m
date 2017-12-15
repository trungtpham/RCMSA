function s = med_scale_estimator(x, th)

min_sig = th*0.5;
max_sig = th*2;

% Remove certain outliers
x(x>max_sig)  = [];

x = abs(x - median(x));
s = 1.4826*mad(x,1);
s = s*2.5;
s = max(min_sig, min(s,max_sig));
end