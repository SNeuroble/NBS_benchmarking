function [y_predict,residuals,x_windowed,y_windowed,y_windowed_std,SSE,residual_r,residual_p]=fit_spline(x,y,smoothing,window_sz)  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function fits a cubic spline for the "Cluster Failure or Power
% Failure" project. Fits a cubic spline to a sliding window avg with 50% 
% overlap.
% Usage: called from summarize_tprs.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sliding window to get smooth model
percent_overlap=0.5; % percent overlap (0.5=50%)
overlap=window_sz*percent_overlap;
nwindows=ceil( (max(x) - min(x) - overlap) / (window_sz - overlap) );
for i=1:nwindows
    x_windowed(i)= min(x)+i*(window_sz-overlap);
    y_in_window=y(x>=(x_windowed(i)-overlap) & x<(x_windowed(i)+overlap));
    y_windowed(i)=mean(y_in_window);
    y_windowed_std(i)=std(y_in_window);
end

% fit spline to smooth model
y_predict = csaps(double(x_windowed),double(y_windowed),smoothing,x);

% descriptive stuff
residuals=y-y_predict;
SSE=sum(residuals.^2);
[residual_r,residual_p]=corr(residuals,x);


