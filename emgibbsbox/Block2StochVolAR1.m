function [h, hshock, h0, kai2States] = Block2StochVolAR1(logy2, N11, ndx11, N22, ndx22, h, rho, hsqrtvcv, KSC, KSCt, Ny, T, rndStream)
% CommonStochVolAR1 performs a Gibbs updating step on a common SV model
% with AR1 SV
%
% Uses Kim, Shephard and Chib normal mixtures
%
%
% See also abcDisturbanceSmoothingSampler1draw, vectorRWsmoothingsampler1draw, getKSC7values, getKSC10values

%   Coded by  Elmar Mertens, em@elmarmertens.com

% note distinction between Ny observables, and Nsv=1 SV processes

% TODO: use sqrthvcv instead of hvcv

% note: assigning N11 and ndx11 directly (to avoid confusion over whether ndx11 could be numeric or logical index)

%% blow up h according to block structure
Nblocks         = 2;
hBlock          = NaN(Ny,T);
hBlock(ndx11,:) = repmat(h(1,:), N11, 1);
hBlock(ndx22,:) = repmat(h(2,:), N22, 1);

%% draw mixture states
% zdraws are standardized draws for each component of the normal mixture
% zdraws is thus Ny x T x Nmixtures
zdraws      = (logy2 - hBlock - KSCt.mean) ./ KSCt.vol;

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

% TODO: implement PS
% vecobs         = obs(:);
% noisevol       = KSC.vol(kai2States(:));
% [h, hshock]    = commonAR1noisePrecisionBasedSampler(vecobs, Ny, T, rho, hvol, noisevol, 1, rndStream);

% DK sampler (quite a bit slower than precision-based sampler):
A     = diag(rho);
B     = hsqrtvcv;
C     = zeros(Ny,Nblocks);
C(ndx11,1) = 1;
C(ndx22,2) = 1;

sqrtR = zeros(Ny,Ny,T);
for n = 1 : Ny
    sqrtR(n,n,:) = KSC.vol(kai2States(n,:));
end

% initial scale levels diffuse, mean fixed at zero
sqrtVh0 = 100 .* eye(Nsv);
Eh0     = zeros(Nsv,1);
[h, hshock, h0] = a2b2c2DisturbanceSmoothingSampler1draw(A, B, C, obs, Eh0, sqrtVh0, ...
        sqrtR, rndStream);


