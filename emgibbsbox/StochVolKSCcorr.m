function [h, h0, hshock, kai2States] = StochVolKSCcorr(logy2, h, hVCV, Eh0, Vh0, KSC, KSCt, Nsv, T, rndStream)
% StochVolKSC performs a Gibbs updating step on a SV model and it works over 
% Nsv independent SV residuals
%
% Uses Kim, Shephard and Chib normal mixtures
%
% USAGE:[h, h0, kai2States] = StochVolKSCcorr(logy2, h, hVCV, Eh0, Vh0, KSC, KSCt, Nsv, T, rndStream)
%
% multivariate case with correlated shocks and RW dynamics
%
% See also rwnoisePrecisionBasedSampler, getKSC7values, getKSC10values

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
% zdraws      = bsxfun(@minus, logy2 - h, KSCt.mean) ./ KSCt.vol;
zdraws      = (logy2 - h - KSCt.mean) ./ KSCt.vol;

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
obs         = logy2 - KSC.mean(kai2States);
 

% precision based sampler
vecobs         = obs(:);
noisevol       = KSC.vol(kai2States(:));
hVCVsqrt       = chol(hVCV, 'lower');
[h, hhat]      = rwnoisePrecisionBasedSampler(vecobs, Nsv, T, hVCVsqrt, noisevol, Eh0, sqrtVh0, 1, rndStream);
        
h0     = hhat(:,1) + hVCVsqrt * randn(rndStream,Nsv,1); % backward simulation
hshock = diff([h0, h], [], 2);




