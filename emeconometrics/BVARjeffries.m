function [PAIhat, sqrtXX, sqrtSSR, dof, R2] = BVARjeffries(Y,X)
% BVARJEFF performs BVAR estimation with jeffries prior (or of dummy obs within Y and X)
% 
% usage: [PAI, sqrtXX, sqrtSSR, dof, R2] = BVARjeffries(Y,X)
%
% Posterior is (PAI, SIGMA) sim MNIW (PAIhat, inv(X'X), Shat, dof)
% with X'X = sqrtXX * sqrtXX'
% and Shat = sqrtSSR * sqrtSSR' / (T - K) where K=size(X,2)
%
%  See also: VARls, dummyobs4BVAR, BVARdraws

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

PAIhat = (sqrtYX / sqrtXX)';
% PAI = (sqrtYX * sqrtInvXX)';

dof = T - K;

% report R2
if nargout > 4
    R2 = 1 - diag(sqrtSSR * sqrtSSR') ./ diag(Y' * Y);
end