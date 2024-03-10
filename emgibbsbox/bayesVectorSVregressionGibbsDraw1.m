function [A, residDraw]  = bayesVectorSVregressionGibbsDraw1(Y, X, iSigmaResid, iVa0a0, iVa0, rndStream)
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

[T, Ny]     = size(Y);
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

iVa = iVa0;
for t = 1 : T
    iVa = iVa + kron(iSigmaResid(:,:,t), X(t,:)' * X(t,:));
end

% posterior variance
sqrt_aSigma    = chol(iVa) \ Ia;            % notice: Matlab's choleski delivers UPPER triangular matrix

% posterior mean
Ytilde = zeros(size(Y));
for t = 1 : T
    Ytilde(t,:) = Y(t,:) * iSigmaResid(:,:,t);
end
XYiSig   = X' * Ytilde;

aTilde   = sqrt_aSigma' * (iVa0a0 + XYiSig(:));

% draw from posterior
aDraw   = sqrt_aSigma * (aTilde + z); % chol_aSigma is the UPPER triangular factorization of aSigma, but this is OK for drawing RV

A       = reshape(aDraw, Nx, Ny);

if nargout > 1
    residDraw   = Y - X * A;
end
