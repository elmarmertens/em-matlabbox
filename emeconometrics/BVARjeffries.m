function [PAI, sqrtXX, sqrtSSR, dof, R2] = BVARjeffries(Y,X)
% BVARJEFF ... 
%  
%   ... 

% get dimensions
[T, K] = size(X);
[~, N] = size(Y);

% perform OLS via QR decomposition
[~,RR]  = qr([X Y],0);
RR      = RR';

% select output
sqrtXX    = RR(1:K,1:K);

sqrtYX    = RR(K+(1:N),1:K);
sqrtSSR   = RR(K+(1:N),K+(1:N));

PAI = (sqrtYX / sqrtXX)';
% PAI = (sqrtYX * sqrtInvXX)';

dof = T - K;

% report R2
if nargout > 4
    R2 = 1 - diag(sqrtSSR * sqrtSSR') ./ diag(Y' * Y);
end