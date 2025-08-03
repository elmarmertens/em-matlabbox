function xbar = ewma(x, lambda, x0)
% EWMA ... 
%  
%   ... 

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 03-Oct-2014 15:07:44 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 8.1.0.604 (R2013a) 
% FILENAME  : ewma.m 


if nargin < 3
    x0 = x(1,:);
end
if isscalar(x0)
    x0 = repmat(x0, 1, size(x,2));
end


xorg = x;
xbar = NaN(size(x));
T    = size(x,1);

t = 1;
nanny = isnan(x(t,:));
if any(nanny) 
    x(t,nanny) = x0(nanny);
end
xbar(t,:) = (1 - lambda) * x(t,:) + lambda * x0;

for t = 2 : T
    nanny = isnan(x(t,:));
    if any(nanny)
        x(t,nanny) = xbar(t-1,nanny);
    end
    xbar(t,:) = (1 - lambda) * x(t,:) + lambda * xbar(t-1,:);
end  

xbar(isnan(xorg)) = NaN; 
