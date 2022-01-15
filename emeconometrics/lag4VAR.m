function [X, Y] = lag4VAR(Y, p, constFlag)
% [X, Y] = lag4VAR(Y, p, constFlag)

%   Coded by  Elmar Mertens, em@elmarmertens.com

% Elmar Mertens
% www.elmarmertens.ch

if nargin < 3
   constFlag = true;
end

[nobs, N] = size(Y);
T = nobs - p;
k = N * p;  

if constFlag
   k = N * p + 1;  
end

% construct regressors
X = ones(T, k);
for i = 1 : p
   X(:, (i - 1) * N + (1 : N))   = Y((p+1 : end) - i,       :);
end
if nargout > 1
   Y = Y(p+1 : end,:);
end
