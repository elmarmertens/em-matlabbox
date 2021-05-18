function [paiDraws, resid, pai, paiSigma]  = bayesRegSysCONSTgibbsDraw(Y, X, Ny, K, ~, SIGMA, iV0pai, iV0, rndStream)
%
%
% See also bayesVectorRegressionGibbsDraw2Steps, bayesVectorRegressionGibbsDraw, lag4VAR

%   Coded by  Elmar Mertens, em@elmarmertens.com


if nargin < 9 || isempty(rndStream)
    rndStream = getDefaultStream;
end

% [T, Ny]  = size(Y);
Nx       = size(X, 2);
Npai      = Ny * K;
Ipai      = eye(Npai);

%% posterior for coefficients
invSIGMA = SIGMA \ eye(Ny);
iV       = iV0 + kron(invSIGMA, X' * X);

cholV    = chol(iV) \ Ipai;  
paiSigma = cholV * cholV';

% posterior mean
pai      = paiSigma * (iV0pai + kron(invSIGMA, X)' * Y(:));

% draw from posterior
paiDraws   = bsxfun(@plus, pai, cholV * randn(rndStream, Npai, 1)); 
paiDraws   = reshape(paiDraws, Nx, Ny);
if nargout > 1
    resid    = Y - X * paiDraws;
end

