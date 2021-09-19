function m = cumnanmedian(x)
% CUMNANMEDIAN ... 
%  
%   ... 

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 02-Feb-2021 22:51:18 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 9.9.0.1570001 (R2020b) Update 4 
% FILENAME  : cumnanmedian.m 


% assume x is vector
x = x(:);

ndx = ~isnan(x);

xtilde = x(ndx);
mtilde = NaN(size(xtilde));
for t = 1 : length(mtilde)
    mtilde(t) = median(xtilde(1:t));
end

m = NaN(size(x));
m(ndx) = mtilde;