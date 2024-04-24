function [A, residDraw, meanA]  = bayesVectorRegressionGibbsDraw1(Y, X, iSigmaResid, iVa0a0, iVa0, rndStream)
% bayesVectorRegressionGibbsDraw performs Gibbs step for vector linear regression model with known variance
%
% [A, residDraw]  = bayesVectorRegressionGibbsDraw1(Y, X, iSigmaResid, iVa0a0, iVa0, rndStream)
%
% iVa0   is inverse prior variance
% iVa0a0 is product of inverse prior variance and prior mean

%   Coded by  Elmar Mertens, em@elmarmertens.com

%% collect input arguments
if nargin < 6 || isempty(rndStream)
    rndStream = getDefaultStream;
end

[~, Ny]     = size(Y);
[~, Nx]     = size(X);

%% draw random numbers as needed
if isnumeric(rndStream)
    z = rndStream;
else
    Na  = length(iVa0a0);
    z   = randn(rndStream, Na, 1);
end


%% posterior for coefficients

% posterior variance
iVa            = iVa0 + kron(iSigmaResid, X' * X);
invsqrtVa      = chol(iVa);   % notice: Matlab's choleski delivers UPPER triangular matrix

% posterior mean
XYiSig   = X' * Y * iSigmaResid;
aTilde   = transpose(invsqrtVa) \ (iVa0a0 + XYiSig(:));
aDraw    = invsqrtVa \ (aTilde + z); 

A       = reshape(aDraw, Nx, Ny);

if nargout > 1
    residDraw   = Y - X * A;
end
if nargout > 2
    meanA = invsqrtVa \ aTilde;
end
