function [A, residDraw, meanA]  = bayesVectorSVregressionGibbsDraw1(Y, X, iSigmaResid, iVa0a0, iVa0, rndStream)
% bayesVectorRegressionGibbsDraw performs Gibbs step for vector linear regression model with known variance
%
% [A, residDraw, meanA]  = bayesVectorSVregressionGibbsDraw1(Y, X, iSigmaResid, iVa0a0, iVa0, rndStream)
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

%% draw random numbers as needed
if isnumeric(rndStream)
    z = rndStream;
else
    Na = length(iVa0a0);
    z  = randn(rndStream, Na, 1);
end


%% posterior for coefficients

% posterior variance
iVa = iVa0;
Xprime = transpose(X); % for better memeory management in the next for loop
for t = 1 : T
    iVa = iVa + kron(iSigmaResid(:,:,t), Xprime(:,t) * transpose(Xprime(:,t)));
end


% posterior mean
Yprime = transpose(Y); % for better memeory management in the next for loop
for t = 1 : T
    Yprime(:,t) = iSigmaResid(:,:,t) * Yprime(:,t);
end
XYiSig   = Xprime * transpose(Yprime);

choliVa  = chol(iVa);
aTilde   = transpose(choliVa) \ (iVa0a0 + XYiSig(:));

% draw from posterior
aDraw   = choliVa \ (aTilde + z); % chol_aSigma is the UPPER triangular factorization of aSigma, but this is OK for drawing RV

A       = reshape(aDraw, Nx, Ny);

if nargout > 1
    residDraw   = Y - X * A;
end
if nargout > 2
    meanA = invsqrtVa \ aTilde;
end
