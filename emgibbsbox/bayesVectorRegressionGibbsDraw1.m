function [A, residDraw]  = bayesVectorRegressionGibbsDraw1(Y, X, iSigmaResid, iVa0a0, iVa0, rndStream)
% bayesVectorRegressionGibbsDraw performs Gibbs step for vector linear regression model with known variance

% iVa0   is inverse prior variance
% iVa0a0 is product of inverse prior variance and prior mean

%   Coded by  Elmar Mertens, em@elmarmertens.com

%% collect input arguments
if nargin < 6 || isempty(rndStream)
    rndStream = getDefaultStream;
end

[~, Ny]     = size(Y);
[~, Nx]     = size(X);
Na          = length(iVa0a0);
Ia          = eye(Na);

%% draw random numbers as needed
if isnumeric(rndStream)
    z = rndStream;
else
    z       = randn(rndStream, Na, 1);
end


%% posterior for coefficients

% posterior variance
iVa            = iVa0 + kron(iSigmaResid, X' * X);
sqrt_aSigma    = chol(iVa) \ Ia;            % notice: Matlab's choleski delivers UPPER triangular matrix

% posterior mean
XYiSig   = X' * Y * iSigmaResid;

aTilde   = sqrt_aSigma' * (iVa0a0 + XYiSig(:));

% draw from posterior
aDraw   = sqrt_aSigma * (aTilde + z); % chol_aSigma is the UPPER triangular factorization of aSigma, but this is OK for drawing RV

A       = reshape(aDraw, Nx, Ny);

if nargout > 1
    residDraw   = Y - X * A;
end
