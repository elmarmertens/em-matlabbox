function [h, hshock, h0, kai2States] = CommonStochVolAR1(logy2, h, rho, hvol, KSC, KSCt, Ny, T, rndStream)
% CommonStochVolAR1 performs a Gibbs updating step on a common SV model
% with AR1 SV
%
% Uses Kim, Shephard and Chib normal mixtures
%
%
% See also abcDisturbanceSmoothingSampler1draw, vectorRWsmoothingsampler1draw, getKSC7values, getKSC10values

%   Coded by  Elmar Mertens, em@elmarmertens.com

% note distinction between Ny observables, adn Nsv=1 SV processes


%% draw mixture states
% zdraws are standardized draws for each component of the normal mixture 
% zdraws is thus Ny x T x Nmixtures
zdraws      = (logy2 - h - KSCt.mean) ./ KSCt.vol;

% construct CDF
% factor of sqrt(2 * pi) can be ommitted for kernel
pdfKernel           = KSCt.pdf ./ KSCt.vol .* exp(-.5 * zdraws.^2); 
cdf                 = cumsum(pdfKernel, 3);                % integrate
cdf(:,:,1:end-1)    = cdf(:,:,1:end-1) ./ cdf(:,:, end); % using automatic expansion 
cdf(:,:,end)        = 1;    % normalize

% draw states
kai2States  = sum(rand(rndStream, Ny, T) > cdf, 3) + 1;


%% KSC State Space
obs   = logy2 - KSC.mean(kai2States);

vecobs          = obs(:);
noisevol        = KSC.vol(kai2States(:));
[h, hshock, h0] = commonAR1noisePrecisionBasedSampler(vecobs, Ny, T, rho, hvol, noisevol, 1, rndStream);

% DK sampler (quite a bit slower than precision-based sampler):
% A     = rho;
% B     = hvol;
% C     = ones(Ny,1);
% sqrtR = zeros(Ny,Ny,T);
% for n = 1 : Ny
%     sqrtR(n,n,:) = KSC.vol(kai2States(n,:));
% end
% sqrtVh0 = 0;
% Eh0     = 0;
% [h, hshock] = a2b2c2DisturbanceSmoothingSampler1draw(A, B, C, obs, Eh0, sqrtVh0, ...
%     sqrtR, rndStream); 



