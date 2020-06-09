function [h, h0, SV, outlierlog2Draws, outlierProb, outlierScaleDraws] = StochVolOutlierKSC(logy2, h, hInno, Eh0, Vh0, outlierlog2Draws, outlierProb, outlieralpha, outlierbeta, outlierStates, KSC, KSCt, Nsv, T, rndStream)
% StochVolOutlierKSC combines KSC Gibbs Sampling for SV with outlier model of Stock-Watson (2016, REStat)
%
% Uses Kim, Shephard and Chib normal mixtures
%
% USAGE: [h, h0, SV, outlierlog2Draws, outlierProb, outlierScaleDraws] = ...
%      StochVolOutlierKSC(logy2, h, hInno, Eh0, Vh0, ...
%      outlierlog2Draws, outlierProb, outlieralpha, outlierbeta, outlierStates, ...
%      KSC, KSCt, Nsv, T, rndStream)
%
%
% See also getKSCvalues

%   Coded by  Elmar Mertens, em@elmarmertens.com



if Nsv > 1
    error('Code supports only Nsv = 1 at this point')
end

if isscalar(Eh0)
    Eh0 = repmat(Eh0, Nsv, 1);
end
if isscalar(Vh0)
    Vh0 = repmat(Vh0, Nsv, 1);
end

%% draw mixture states
% zdraws are standardized draws for each component of the normal mixture 
% zdraws is thus Nsv x T x 7 
zdraws      = bsxfun(@minus, logy2 - h - outlierlog2Draws, KSCt.mean) ./ KSCt.vol;

% construct CDF
% factor of sqrt(2 * pi) can be ommitted for kernel

pdfKernel           = KSCt.pdf ./ KSCt.vol .* exp(-.5 * zdraws.^2); 
cdf                 = cumsum(pdfKernel, 3);                % integrate
cdf(:,:,1:end-1)    = bsxfun(@rdivide, cdf(:,:,1:end-1), cdf(:,:, end)); 
cdf(:,:,end)        = 1;    % normalize

% draw states
kai2States  = sum(bsxfun(@gt, rand(rndStream, Nsv, T), cdf), 3) + 1;
obs         = logy2 - KSC.mean(kai2States) - outlierlog2Draws;

%% KSC State Space
h = NaN(Nsv,T);
h0 = NaN(Nsv,1);

for n = 1 : Nsv
    [h(n,:), h0(n)] = smoothingsamplerRWnoise(obs(n,:),hInno(n)^2,KSC.var(kai2States(n,:)),...
        Eh0(n),Vh0(n),rndStream);
end

%% outlier PDF
% outlierPdf is Nsurvey times T  times Nstates
% outlierPdf2 = cat(3, repmat(1 - outlierProb, 1, T), bsxfun(@times, outlierProb, repmat(1 / outlierNgrid, 1, T, outlierNgrid)));
outlierPdf  = cat(3, repmat(1 - outlierProb, 1, T), repmat(outlierProb / outlierStates.Ngrid, 1, T, outlierStates.Ngrid));
% checkdiff(outlierPdf, outlierPdf2);
% keyboard

%% outlier states
edraws      = bsxfun(@minus, logy2 - h - KSC.mean(kai2States), permute(outlierStates.log2values, [1 3 2]));
zdraws      = bsxfun(@rdivide, edraws, KSC.vol(kai2States));

% pdfKernel   = exp(-.5 * zdraws.^2);  
% division by KSC.vol is unnecessarty for this kernel, since same vol would apply across outlierStates
% pdfKernel   = bsxfun(@rdivide, pdfKernel, KSC.vol(kai2States));

% check: bsxfun needed?
% pdfKernel   = bsxfun(@times, outlierPdf, pdfKernel);
pdfKernel   = outlierPdf .* exp(-.5 * zdraws.^2);

cdf                 = cumsum(pdfKernel, 3);                % integrate
cdf(:,:,1:end-1)    = bsxfun(@rdivide, cdf(:,:,1:end-1), cdf(:,:,end)); 
cdf(:,:,end)        = 1;    % normalize


% draw states
ndx               = sum(bsxfun(@gt, rand(rndStream, Nsv, T), cdf), 3) + 1;
outlierlog2Draws  = outlierStates.log2values(ndx);
outlierScaleDraws = outlierStates.values(ndx);

%% update outlierProb
Noutlier    = sum(ndx > 1, 2);
alpha       = outlieralpha + Noutlier;
beta        = outlierbeta + (T - Noutlier);
for n = 1 : Nsv
    outlierProb(n) = betadraw(alpha(n), beta(n), 1, rndStream);
    %     outlierProb(n) = betarnd(alpha(n), beta(n), 1);
end

%% construct SV
SV = exp((h + outlierlog2Draws) / 2);
