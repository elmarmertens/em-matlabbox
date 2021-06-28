function [h, hbar, hresid, htilde, kai2States] = StochVolKSCcorrAR1(logy2, h, rho, hVCV, Eh0, Vh0, KSC, KSCt, Nsv, T, rndStream)
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
obs         = logy2 - KSC.mean(kai2States);

%% AR1 parameters for SV (fixed)
% rho = repmat(0.99, Nsv, 1);

%% KSC State Space
zerosNsv  = zeros(Nsv);
Insv      = eye(Nsv);
A     = [diag(rho) zerosNsv; zerosNsv Insv];
B     = [chol(hVCV)'; zerosNsv];
C     = [Insv Insv];
sqrtR = zeros(Nsv,Nsv,T);
for n = 1 : Nsv
    sqrtR(n,n,:) = KSC.vol(kai2States(n,:));
end

sqrtVh0 = diag(sqrt(Vh0));

sqrtVhtilde  = zeros(Nsv); % Note: need fixed prior, not depended on estimated rhos (alt: use prior rho)
x0           = [zeros(Nsv, 1); Eh0];
sqrtVx0      = [sqrtVhtilde, zerosNsv; zerosNsv sqrtVh0];
[H, Hshock, H0] = abcDisturbanceSmoothingSampler1draw(A, B, C, obs, x0, sqrtVx0, ...
    sqrtR, rndStream); 

h      = H(1:Nsv,:) + H(Nsv+1:end,:); % C * H
hbar   = H0(Nsv+1:end);
htilde = cat(2, H0(1:Nsv,:), H(1:Nsv,:));
hresid = Hshock(1:Nsv,:);
