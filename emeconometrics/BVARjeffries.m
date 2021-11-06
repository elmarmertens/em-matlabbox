function [PAIhat, sqrtXX, sqrtSSR, dof, postllf, R2] = BVARjeffries(Y,X)
% BVARJEFF performs BVAR estimation with jeffries prior (or of dummy obs within Y and X)
% 
% usage: [PAI, sqrtXX, sqrtSSR, dof, postllf, R2] = BVARjeffries(Y,X)
%
% Posterior is (PAI, SIGMA) sim MNIW (PAIhat, inv(X'X), Shat, dof)
% with X'X = sqrtXX * sqrtXX'
% and Shat = sqrtSSR * sqrtSSR' / (T - K) where K=size(X,2)
%
% postllf:
% - is contribution of posterior to log-like 
% - to compute llf call [~,~,~,~,priorllf] = BVARjeffries(Ystar, Xstar) where Ystar and Xstar are dummy obs for prior
% - note: requires proper prior with Tstar > K
% - llf = postlf - priorllf
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

% post-loof
if nargout > 4
    logdetS  = 2 * sum(log(abs(diag(sqrtSSR))));
    logdetXX = 2 * sum(log(abs(diag(sqrtXX))));
    postllf = - N * T / 2 * log(2 * pi) -dof/2 * logdetS - N/2 * logdetXX + N * dof /2 + sum(gammaln((dof + 1 - 1:N) / 2));

    checkdiff(logdetS, log(det(sqrtSSR * sqrtSSR')));
    checkdiff(logdetXX, log(det(X' * X)));
end
% report R2
if nargout > 5
    R2 = 1 - diag(sqrtSSR * sqrtSSR') ./ diag(Y' * Y);
end