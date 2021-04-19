function [aDraws, resid, a, aSigma]  = bayesRegSysSVgibbsDraw(Y, X, iSigmaResid, iVa0, iV0, rndStream)
%
%
% See also bayesVectorRegressionGibbsDraw2Steps, bayesVectorRegressionGibbsDraw, lag4VAR

%   Coded by  Elmar Mertens, em@elmarmertens.com

%% VERSION INFO
% AUTHOR    : Elmar Mertens
% $DATE     : 01-Sep-2010 11:35:11 $
% $Revision : 1.00 $
% DEVELOPED : 7.8.0.347 (R2009a)
% FILENAME  : bayesVARgibbsDraw.m

if nargin < 6 || isempty(rndStream)
    rndStream = getDefaultStream;
end

[T, Ny]  = size(Y);
Nx       = size(X, 2);
Na       = length(iVa0);
Ia       = eye(Na);

%% posterior for coefficients

% posterior variance
inv_aSigma     = iV0;
for t = 1 : T
    inv_aSigma = inv_aSigma + kron(iSigmaResid(:,:,t), X(t,:)' * X(t,:));
end
chol_aSigma    = chol(inv_aSigma) \ Ia;            % notice: Matlab's choleski delivers UPPER triangular matrix
aSigma         = chol_aSigma * chol_aSigma';       % checkdiff(aSigma, inv(inv_aSigma));

% posterior mean
Ytilde = zeros(size(Y));
for t = 1 : T
    Ytilde(t,:) = Y(t,:) * iSigmaResid(:,:,t);
end
tmp      = X' * Ytilde;
a        = aSigma * (iVa0 + tmp(:));

% draw from posterior
aDraws   = bsxfun(@plus, a, chol_aSigma * randn(rndStream, Na, 1)); % chol_aSigma is the UPPER triangular factorization of aSigma, but this is OK for drawing RV
aDraws   = reshape(aDraws, Nx, Ny);
if nargout > 1
    resid    = Y - X * aDraws;
end

