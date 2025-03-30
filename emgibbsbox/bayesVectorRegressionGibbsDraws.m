function [PAI, residDraw]  = bayesVectorRegressionGibbsDraws(Y, X, iSigmaResid, iVpai0pai0, iVpai0, Ndraws, rndStream)
% bayesVectorRegressionGibbsDraw performs Gibbs step for vector linear regression model with known variance
%
% [PAI, residDraw]  = bayesVectorRegressionGibbsDraws(Y, X, iSigmaResid, iVpai0pai0, iVpai0, Ndraws, rndStream)
%
% iVpai0   is inverse prior variance
% iVpai0pai0 is product of inverse prior variance times prior mean
% rndStream can be stream or a set of standard-normal random numbers

%   Coded by  Elmar Mertens, em@elmarmertens.com

%% collect input arguments
if nargin < 6 || isempty(Ndraws)
    Ndraws = 1;
end
if nargin < 7 || isempty(rndStream)
    rndStream = getDefaultStream;
end

Ny     = size(Y,2);
Nx     = size(X,2);
Npai   = length(iVpai0pai0);
Ipai   = eye(Npai);

%% draw random numbers as needed
if isnumeric(rndStream)
    z = rndStream;
else
    z = randn(rndStream, Npai, Ndraws);
end


%% posterior for coefficients

% posterior variance
iVpai          = iVpai0 + kron(iSigmaResid, X' * X);
sqrt_paiSigma  = chol(iVpai) \ Ipai;            % notice: Matlab's choleski delivers UPPER triangular matrix

% posterior mean
XYiSig   = X' * Y * iSigmaResid;

paiTilde  = sqrt_paiSigma' * (iVpai0pai0 + XYiSig(:));
paiDraw   = sqrt_paiSigma * (paiTilde + z); % chol_aSigma is the UPPER triangular factorization of aSigma, but this is OK for drawing RV

PAI       = reshape(paiDraw, Nx, Ny, Ndraws);

if nargout > 1
    % residDraw   = NaN(T,Ny,Ndraws);
    % for n = 1 : Ndraws % todo: use pagemetimes
    %     residDraw(:,:,n) = Y - X * PAI(:,:,n);
    % end

    residDraw = Y - pagemtimes(X, PAI);

end
