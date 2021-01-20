function r = reless(w)
% RELESS ... 
%  
%   ... 

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 16-Jan-2021 22:11:58 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 9.9.0.1538559 (R2020b) Update 3 
% FILENAME  : reless.m 

w = w(:);
N = size(w,1);
r = 1 ./ sum(w.^2,1) / N;
