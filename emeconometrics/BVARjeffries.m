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
% - log marginal data density of the supplied sample under the Jeffreys
%   prior p(PAI,SIGMA) propto |SIGMA|^(-(N+1)/2):
%       log m(Y) = -N*dof/2*log(pi) - N/2*log|X'X| - dof/2*log|S|
%                  + log GammaN(dof/2),   dof = T - K
%   with GammaN the multivariate gamma function
% - to compute llf call [~,~,~,~,priorllf] = BVARjeffries(Ystar, Xstar) where Ystar and Xstar are dummy obs for prior
% - note: requires proper prior with Tstar >= K + N (dof >= N)
% - llf = postllf - priorllf
%
%  See also: VARls, dummyobs4BVAR, BVARdraws, unittestBVARjeffries


% get dimensions
[T, K] = size(X);
[~, N] = size(Y);

if T < K + N
    error('em:BVARjeffries', 'need T >= K + N (T=%d, K=%d, N=%d)', T, K, N)
end

% perform OLS via QR decomposition
RR      = qr([X Y],'econ'); % single output skips forming Q
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
    logGammaN = N*(N-1)/4 * log(pi) + sum(gammaln((dof + 1 - (1:N)) / 2));
    postllf  = - N*dof/2 * log(pi) - dof/2 * logdetS - N/2 * logdetXX + logGammaN;
end
% report R2
if nargout > 5
    R2 = 1 - sum(sqrtSSR.^2, 2) ./ sum(Y.^2, 1)';
end

end % function BVARjeffries

