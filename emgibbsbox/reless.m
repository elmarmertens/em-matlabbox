function r = reless(w)
% RELESS computes relative effective sample size
%
% USAGE r = reless(w)
%
% w can be N-element vector or a T x N matrix of weights (T observations of N weights)
%  


%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 16-Jan-2021 22:11:58 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 9.9.0.1538559 (R2020b) Update 3 
% FILENAME  : reless.m 

if scalar(w)
    w = w(:);
end
N = size(w,1);
r = 1 ./ sum(w.^2,1) / N;
