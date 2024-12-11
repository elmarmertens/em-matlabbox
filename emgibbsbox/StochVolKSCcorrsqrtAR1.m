function [h, hbar, hshock, htilde, kai2States] =  ...
     StochVolKSCcorrsqrtAR1(logy2, h, rho, hVCVsqrt, Eh0, sqrtVh0, KSC, KSCt, Nsv, T, rndStream)
% StochVolKSCcorrsqrtAR1 performs a Gibbs updating step on a SV model with AR1
% dynamics and correlated shocks
%
% Uses Kim, Shephard and Chib normal mixtures
%
%
% See also abcDisturbanceSmoothingSampler1draw, vectorRWsmoothingsampler1draw, getKSC7values, getKSC10values

%   Coded by  Elmar Mertens, em@elmarmertens.com

if isscalar(Eh0)
    Eh0 = repmat(Eh0, Nsv, 1);
end
if isscalar(sqrtVh0)
    sqrtVh0 = sqrtVh0 * speye(Nsv);
end


%% draw mixture states
% zdraws are standardized draws for each component of the normal mixture 
% zdraws is thus Nsv x T x Nmixtures
zdraws      = (logy2 - h - KSCt.mean) ./ KSCt.vol;

% construct CDF
% factor of sqrt(2 * pi) can be ommitted for kernel
pdfKernel           = KSCt.pdf ./ KSCt.vol .* exp(-.5 * zdraws.^2); 
cdf                 = cumsum(pdfKernel, 3);                % integrate
cdf(:,:,1:end-1)    = cdf(:,:,1:end-1) ./ cdf(:,:, end); % using automatic expansion 
cdf(:,:,end)        = 1;    % normalize

% draw states
kai2States  = sum(rand(rndStream, Nsv, T) > cdf, 3) + 1;


%% KSC State Space
obs       = logy2 - KSC.mean(kai2States);
noisevol  = KSC.vol(kai2States);

[h, hbar, hshock, htilde] = sampleVAR1noisePBS(obs, Nsv, T, rho, hVCVsqrt, Eh0, sqrtVh0, noisevol, rndStream);
