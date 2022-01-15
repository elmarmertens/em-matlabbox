function [h, hbar, htilde, hresid, kai2States] = ...
    StochVolAR1(logy2, h, rho, hVCVsqrt, Eh0, sqrtVh0, ...
	KSC, KSCt, Nsv, T, rndStream)
% StochVolKSCcorrAR1 performs a Gibbs updating step on a SV model with AR1
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
    sqrtVh0 = sqrtVh0 * eye(Nsv);
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
obs       = logy2 - KSC.mean(kai2States);
zerosNsv  = zeros(Nsv);
Insv      = eye(Nsv);
A     = [diag(rho) zerosNsv; zerosNsv Insv];
B     = [hVCVsqrt; zerosNsv];
C     = [Insv Insv];
sqrtR = zeros(Nsv,Nsv,T);
for n = 1 : Nsv
    sqrtR(n,n,:) = KSC.vol(kai2States(n,:));
end

sqrtVhtilde  = zeros(Nsv); % Note: need fixed prior, not depended on estimated rhos (alt: use prior rho)
x0           = [zeros(Nsv, 1); Eh0];
sqrtVx0      = [sqrtVhtilde, zerosNsv; zerosNsv sqrtVh0];
[H, Hshock, H0] = a2b2c2DisturbanceSmoothingSampler1draw(A, B, C, obs, x0, sqrtVx0, ...
    sqrtR, rndStream); 

h      = H(1:Nsv,:) + H(Nsv+1:end,:); % C * H
hbar   = H0(Nsv+1:end);
htilde = cat(2, H0(1:Nsv,:), H(1:Nsv,:));
hresid = Hshock(1:Nsv,:);
