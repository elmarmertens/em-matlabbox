function xbar = ewmaND(x, lambda, x0)
% EWMA ... 
%  
%   ... 


%% fold into 2D
xsize = size(x);

x     = reshape(x, xsize(1), []);


%% defaults

if nargin < 3
    x0 = x(1,:);
end
if isscalar(x0)
    x0 = repmat(x0, 1, size(x,2));
end

%% compute ewma
T    = size(x,1);

t = 1;
xbar(t,:) = (1 - lambda) * x(t,:) + lambda * x0;

for t = 2 : T
    xbar(t,:) = (1 - lambda) * x(t,:) + lambda * xbar(t-1,:);
end  

xbar = reshape(xbar, xsize);
