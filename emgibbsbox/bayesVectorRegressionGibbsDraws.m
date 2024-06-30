function [A, residDraw]  = bayesVectorRegressionGibbsDraws(Y, X, iSigmaResid, iVa0a0, iVa0, Ndraws, rndStream)
% bayesVectorRegressionGibbsDraw performs Gibbs step for vector linear regression model with known variance

% iVa0   is inverse prior variance
% iVa0a0 is product of inverse prior variance and prior mean

%   Coded by  Elmar Mertens, em@elmarmertens.com

%% collect input arguments
if nargin < 6 || isempty(Ndraws)
    Ndraws = 1;
end
if nargin < 7 || isempty(rndStream)
    rndStream = getDefaultStream;
end

[T, Ny]     = size(Y);
[~, Nx]     = size(X);
Na          = length(iVa0a0);
Ia          = eye(Na);

%% draw random numbers as needed
if isnumeric(rndStream)
    z = rndStream;
else
    z = randn(rndStream, Na, Ndraws);
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

% check:
% aSigma  = sqrt_aSigma * sqrt_aSigma';       % checkdiff(aSigma, inv(inv_aSigma));
% a       = aSigma * (iVa0a0 + XYiSig(:));
% checkdiff(aDraw, a + sqrt_aSigma * z);

A       = reshape(aDraw, Nx, Ny, Ndraws);

if nargout > 1
    residDraw   = NaN(T,Ny,Ndraws);
    for n = 1 : Ndraws % todo: use pagemetimes
        residDraw(:,:,n) = Y - X * A(:,:,n);
    end
end
