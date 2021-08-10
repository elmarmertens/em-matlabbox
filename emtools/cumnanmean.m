function m = cumnanmean(x)
% CUMNANMEAN ... 
%  
%   ... 

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 02-Feb-2021 22:47:47 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 9.9.0.1570001 (R2020b) Update 4 
% FILENAME  : cumnanmean.m 

% assume x is vector
x = x(:);

ndx = ~isnan(x);

xtilde = x(ndx);
mtilde = cumsum(xtilde) ./ (1:length(xtilde))';

m = NaN(size(x));
m(ndx) = mtilde;