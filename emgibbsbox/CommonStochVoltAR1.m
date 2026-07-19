function [h, hshock, h0, kai2States, tscalelog2Draws, tdofDraw] = CommonStochVoltAR1(y2, logy2, h, rho, hvol, KSC, KSCt, tdof, Ny, T, rndStream)
% CommonStochVolAR1 performs a Gibbs updating step on a common SV model with multivariate t
% with AR1 SV
%
% Uses Kim, Shephard and Chib normal mixtures
%
%
% See also abcDisturbanceSmoothingSampler1draw, vectorRWsmoothingsampler1draw, getKSC7values, getKSC10values, CommonStochVolAR1

%   Coded by  Elmar Mertens, em@elmarmertens.com

% note distinction between Ny observables, and Nsv=1 SV processes

%% t-shocks via outlier draws (following JPR04)
SSR        = sum(y2 .* exp(-h),1); % common scaleSV, thus summing over Ny
SSRgrid    = repmat(permute(SSR, [1 3 2]), [1 tdof.Ndof 1]);


dataloglike   = - .5 * (tdof.values + Ny) .* sum(log(1 + (SSRgrid ./ tdof.values)), 3); 
loglike       = tdof.loglike0 + dataloglike;

% loglikeCheck   = sum(gammaln(.5 * (tdof.values + Ny)) - gammaln(.5 * tdof.values) - .5 * Ny .* log(tdof.values) - .5 * Ny * log(pi) - .5 * (tdof.values + Ny) .* log(1 + (SSRgrid ./ tdof.values)), 3);
% checkdiff(loglike, loglikeCheck);
% 
% loglike0Check   = T * (gammaln(.5 * (tdof.values + Ny)) - gammaln(.5 * tdof.values) - .5 * Ny .* log(tdof.values)  - .5 * Ny * log(pi));
% checkdiff(tdof.loglike0, loglike0Check);


logposteriorKernel = tdof.logprior + loglike;

% subtract const to avoid overflow
logposteriorKernelstar  = logposteriorKernel  - max(logposteriorKernel, [], 2);


cdf            = cumsum(exp(logposteriorKernelstar), 2);
cdf(:,1:end-1) = cdf(:,1:end-1) ./ cdf(:,end);
cdf(:,end)     = 1;
dofState       = sum(rand(rndStream, 1, 1) > cdf, 2) + 1;
tdofDraw       = tdof.values(dofState); % (scalar)


scalePosterior     = tdofDraw + SSR; 
chi2draws          = chi2rnd(tdofDraw + Ny, 1, T); 
tscalelog2Draws    = log(scalePosterior) - log(chi2draws);
% slightly slower: tscalelog2Draws    = log(igamdraw(scalePosterior, tdofDraw + Ny));

%% draw mixture states
% zdraws are standardized draws for each component of the normal mixture 
% zdraws is thus Ny x T x Nmixtures
zdraws      = (logy2 - h - tscalelog2Draws - KSCt.mean) ./ KSCt.vol;

% construct CDF
% factor of sqrt(2 * pi) can be ommitted for kernel
pdfKernel           = KSCt.pdf ./ KSCt.vol .* exp(-.5 * zdraws.^2); 
cdf                 = cumsum(pdfKernel, 3);                % integrate
cdf(:,:,1:end-1)    = cdf(:,:,1:end-1) ./ cdf(:,:, end); % using automatic expansion 
cdf(:,:,end)        = 1;    % normalize

% draw states
kai2States  = sum(rand(rndStream, Ny, T) > cdf, 3) + 1;


%% KSC State Space
obs   = logy2 - KSC.mean(kai2States) - tscalelog2Draws;

vecobs          = obs(:);
noisevol        = KSC.vol(kai2States(:));
[h, hshock, h0] = commonAR1noisePrecisionBasedSampler(vecobs, Ny, T, rho, hvol, noisevol, 1, rndStream);



