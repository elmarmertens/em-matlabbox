function [iwishcholDraw, cholSigmaT, dof] = bayesSQRTVCVgibbsDraw1(Sigma0T, dof0, resid, rndStream)
% bayesSQRTVCVgibbsDraw1 iwishcholDraws SQRT factor of variance-covariance matrix
% [iwishcholDraws, cholSigmaT, dof] = bayesSQRTVCVgibbsDraw1(Sigma0T, dof0, resid, rndStream)
% Prior is (Sigma0T, dof0) invWishart and posterior is (SigmaT, dof) invWishart
% Recall: Mean of inverse Wishart is SigmaT / (dof - N - 1)
% dof = dof0 + T
%
% resid is T x N
%
%
% See also iwishdraw, bayesVCVgibbsDraw1, bayesVCVgibbsDraw

%   Coded by  Elmar Mertens, em@elmarmertens.com


%% parse inputs
if nargin < 4 || isempty(rndStream)
    rndStream = getDefaultStream;
end


if ~isscalar(dof0)
    error('dof must be a scalar')
end

%% draw random normals
if isnumeric(rndStream)
    z = rndStream;
else
    [T, Ny] = size(resid);
    dof     = dof0 + T;
    z       = randn(rndStream, Ny, dof);
end

%% Posterior Update
% Sigma0T     = cholSigma0T * cholSigma0T';
SigmaT      = Sigma0T + resid' * resid;
cholSigmaT  = chol(SigmaT)';

% alt QR code: (way slower)
% [~,sqrtSigmaT] = qr([cholSigma0T'; resid], 0);
% sqrtSigmaT     = sqrtSigmaT';
% checkdiff(SigmaT, sqrtSigmaT * sqrtSigmaT');


%% compute invWishart iwishDraws from iW(Sigma, dof)
sqrtZZ         = chol(z * z'); % note absence of transpose
iwishcholDraw  = cholSigmaT / sqrtZZ;

% wishdraw = cholSigmaT / (z * z') * cholSigmaT';
% checkdiff(wishdraw, iwishcholDraw * iwishcholDraw');
