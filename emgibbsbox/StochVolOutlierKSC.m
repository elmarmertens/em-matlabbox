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
% See also getKSC7values, getKSC10values, StochVolOutlierKSCcorrsqrt

%   Coded by  Elmar Mertens, em@elmarmertens.com


if isscalar(Eh0)
    Eh0 = repmat(Eh0, Nsv, 1);
end
if isscalar(Vh0)
    sqrtVh0 = sqrt(Vh0) * speye(Nsv);
end
if isvector(Vh0)
    sqrtVh0 = sparse(diag(sqrt(Vh0))); % better to define as speye in
                                     % callin function, this is just a backstop
elseif ~issparse(Vh0)
    Vh0 = sparse(Vh0); % better to define as sparse in
                               % calling function, this is just a backstop
    sqrtVh0 = chol(Vh0, 'lower');
end


%% draw mixture states
% zdraws are standardized draws for each component of the normal mixture 
% zdraws is thus Nsv x T x Nmixtures
% zdraws      = bsxfun(@minus, logy2 - h - outlierlog2Draws, KSCt.mean) ./ KSCt.vol;
zdraws      = (logy2 - h - outlierlog2Draws - KSCt.mean) ./ KSCt.vol;

% construct CDF
% factor of sqrt(2 * pi) can be ommitted for kernel
pdfKernel           = KSCt.pdf ./ KSCt.vol .* exp(-.5 * zdraws.^2); 
cdf                 = cumsum(pdfKernel, 3);                % integrate
% cdf(:,:,1:end-1)    = bsxfun(@rdivide, cdf(:,:,1:end-1), cdf(:,:, end)); 
cdf(:,:,1:end-1)    = cdf(:,:,1:end-1) ./ cdf(:,:, end); % using automatic expansion 
cdf(:,:,end)        = 1;    % normalize

% draw states
% kai2States  = sum(bsxfun(@gt, rand(rndStream, Nsv, T), cdf), 3) + 1;
kai2States  = sum(rand(rndStream, Nsv, T) > cdf, 3) + 1;


%% KSC State Space
obs         = logy2 - KSC.mean(kai2States) - outlierlog2Draws;

% precision based sampler
vecobs         = obs(:);
noisevol       = KSC.vol(kai2States(:));
hVCVsqrt       = sparse(diag(hInno));
[h, hhat]      = rwnoisePrecisionBasedSampler(vecobs, Nsv, T, hVCVsqrt, noisevol, Eh0, sqrtVh0, 1, rndStream);
        
h0     = hhat(:,1) + hVCVsqrt * randn(rndStream,Nsv,1); % backward simulation
% hshock = diff([h0, h], [], 2);


% h = NaN(Nsv,T);
% h0 = NaN(Nsv,1);
% 
% for n = 1 : Nsv
%     [h(n,:), h0(n)] = smoothingsamplerRWnoise(obs(n,:),hInno(n)^2,KSC.var(kai2States(n,:)),...
%         Eh0(n),Vh0(n),rndStream);
% end

%% outlier PDF
% outlierPdf is Nsurvey times T  times Nstates
% outlierPdf2 = cat(3, repmat(1 - outlierProb, 1, T), bsxfun(@times, outlierProb, repmat(1 / outlierNgrid, 1, T, outlierNgrid)));
outlierPdf  = cat(3, repmat(1 - outlierProb, 1, T), repmat(outlierProb / outlierStates.Ngrid, 1, T, outlierStates.Ngrid));

%% outlier states
edraws      = bsxfun(@minus, logy2 - h - KSC.mean(kai2States), permute(outlierStates.log2values, [1 3 2]));
zdraws      = bsxfun(@rdivide, edraws, KSC.vol(kai2States));

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
% for n = 1 : Nsv
%     outlierProb(n) = betadraw(alpha(n), beta(n), 1, rndStream);
    % re matlab's betarnd:
    % - does not seem to support randomStreams (but compatible with parfor
    % according to documentation
    % - appears to be faster now (was different in earlier versions)
% end

%% construct SV
SV = exp((h + outlierlog2Draws) / 2);
