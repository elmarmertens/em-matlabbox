function [h, h0, kai2States] = StochVolKSC(logy2, h, hInno, Eh0, Vh0, KSC, KSCt, Nsv, T, rndStream)
% StochVolKSC performs a Gibbs updating step on a SV model and it works over 
% Nsv independent SV residuals
%
% Uses Kim, Shephard and Chib normal mixtures
%
% USAGE: [h, kai2States] =  StochVolKSC(logy2, h, hInno, Eh0, Vh0, KSC, KSCt, Nsv, T, rndStream)
%
% where h are the log-SV's and hInno is the *volatility* in the innovations of h
%
% assumes scalar SV process, see StochVolKSCcorr for multivariate case
%
% See also smoothingsamplerRWnoise, StochVolKSCcorr, getKSC7values, getKSC10values

%   Coded by  Elmar Mertens, em@elmarmertens.com


%% VERSION INFO
% AUTHOR    : Elmar Mertens
% $DATE     : 28-Aug-2009 12:07:01 $
% $Revision : 1.00 $
% DEVELOPED : 7.7.0.471 (R2008b)
% FILENAME  : StochVolKSC.m.m

if isscalar(Eh0)
    Eh0 = repmat(Eh0, Nsv, 1);
end
if isscalar(Vh0)
    Vh0 = repmat(Vh0, Nsv, 1);
end

%% CORRIGENDUM CHANGES ORDER OF GIBBS STEPS!

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
obs = logy2 - KSC.mean(kai2States);
h0 = NaN(Nsv,1);
for n = 1 : Nsv
   [h(n,:), h0(n)] = smoothingsamplerRWnoise(obs(n,:),hInno(n)^2,KSC.var(kai2States(n,:)),...
       Eh0(n),Vh0(n),rndStream);
end

% alt code: call vectorRWsmoothingsampler1draw

