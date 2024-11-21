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
zerosNsv  = zeros(Nsv);
Insv      = eye(Nsv);
A     = [diag(rho) zerosNsv; zerosNsv Insv];
B     = [hVCVsqrt; zerosNsv];
C     = [Insv Insv];
sqrtR = zeros(Nsv,Nsv,T);
for n = 1 : Nsv
    sqrtR(n,n,:) = KSC.vol(kai2States(n,:));
end


x0           = [zeros(Nsv, 1); Eh0];
sqrtVhtilde  = eye(Nsv); % Note: need fixed prior, not dependent on estimated rhos (alt: use prior rho)
sqrtVx0      = [sqrtVhtilde, zerosNsv; zerosNsv sqrtVh0];
[H, Hshock, H0] = a2b2c2DisturbanceSmoothingSampler1draw(A, B, C, obs, x0, sqrtVx0, ...
    sqrtR, rndStream); 

h      = H(1:Nsv,:) + H(Nsv+1:end,:); % C * H
hbar   = H0(Nsv+1:end);
htilde = cat(2, H0(1:Nsv,:), H(1:Nsv,:)); % demeaned component, including lagged value
hshock = Hshock(1:Nsv,:);
