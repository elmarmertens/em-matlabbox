function [paiDraws, resid, pai, paiSigma]  = bayesRegSysSVgibbsDraw(Y, X, Ny, K, T, invA, SV, iV0pai, iV0, rndStream)
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

if nargin < 10 || isempty(rndStream)
    rndStream = getDefaultStream;
end

% [T, Ny]  = size(Y);
Nx       = size(X, 2);
% K        = length(iV0pai);
Npai      = Ny * K;

SVprime   = SV'; % for better memory alignment
Xprime    = X';
Yprime    = Y';

%% posterior for coefficients

inv_paiSigma    = iV0;
Ytilde          = zeros(size(Yprime));
Aprime          = (eye(Ny) / invA)';

for t = 1 : T

    % prepare posterior variance
    sqrtinvSigmaResid = Aprime * diag( 1 ./ SVprime(:,t));
    C                 = kron(sqrtinvSigmaResid, Xprime(:,t));
    inv_paiSigma      = inv_paiSigma + C * C';
    

    % prepare posterior mean
    iSigmaResid = sqrtinvSigmaResid * sqrtinvSigmaResid';
    Ytilde(:,t) = iSigmaResid * Yprime(:,t);

end

Ipai           = eye(Npai);
chol_paiSigma  = chol(inv_paiSigma) \ Ipai;        % notice: Matlab's choleski delivers UPPER triangular matrix
paiSigma       = chol_paiSigma * chol_paiSigma';       % checkdiff(aSigma, inv(inv_aSigma));

% posterior mean
XYtilde  = (Ytilde * X)';
pai      = paiSigma * (iV0pai + XYtilde(:));

% draw from posterior
paiDraws   = bsxfun(@plus, pai, chol_paiSigma * randn(rndStream, Npai, 1)); % chol_aSigma is the UPPER triangular factorization of aSigma, but this is OK for drawing RV
paiDraws   = reshape(paiDraws, Nx, Ny);
if nargout > 1
    resid    = Y - X * paiDraws;
end

