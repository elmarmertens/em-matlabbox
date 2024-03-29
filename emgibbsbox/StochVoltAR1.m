function [h, hbar, htilde, hresid, SV, ...
    SVtScalelog2Draws, SVtdofDraws] = ...
    StochVoltAR1(y2, logy2, h, rho, hVCVsqrt, Eh0, sqrtVh0, ...
    tdof, ...
    KSC, KSCt, Nsv, T, rndStream)
% 
%
%
% See also getKSC7values, getKSC10values, StochVolOutlierKSCcorrsqrt

%   Coded by  Elmar Mertens, em@elmarmertens.com


if isscalar(Eh0)
    Eh0 = repmat(Eh0, Nsv, 1);
end
if isscalar(sqrtVh0)
    sqrtVh0 = sqrtVh0 * eye(Nsv);
end

%% t-shocks via outlier draws (following JPR04)
y2scaled        = y2 .* exp(-h);
y2scaledGrid    = repmat(permute(y2scaled, [1 3 2]), [1 tdof.Ndof 1]);


dataloglike   = - .5 * (tdof.values + 1) .* sum(log(tdof.values + y2scaledGrid), 3); 
loglike       = tdof.loglike0 + dataloglike;

% check against stats pdf
% dofGrid = repmat(tdof.values, [Nsv, 1, T]);
% loglike2 = sum(log(tpdf(sqrt(y2scaledGrid), dofGrid)), 3);
% checkdiff(loglike, loglike2);

logposteriorKernel = tdof.logprior + loglike;
% note: adding prior could be dropped from kernel as long as the prior is uniform
% subtract const to avoid overflow
logposteriorKernelstar  = logposteriorKernel  - max(logposteriorKernel, [], 2);


cdf            = cumsum(exp(logposteriorKernelstar), 2);
cdf(:,1:end-1) = cdf(:,1:end-1) ./ cdf(:,end);
cdf(:,end)     = 1;
dofStates      = sum(rand(rndStream, Nsv, 1) > cdf, 2) + 1;
SVtdofDraws    = tdof.values(dofStates)'; % transpose!


scalePosterior     = SVtdofDraws + y2scaled; % note the explicit vector expansion of SVtdof
% note: matlab doc says stats box handles parallel streams automatically via the global stream ....
chi2draws          = chi2rnd(repmat(SVtdofDraws + 1, 1, T));
SVtScalelog2Draws  = log(scalePosterior) - log(chi2draws);


%% draw mixture states
% zdraws are standardized draws for each component of the normal mixture 
% zdraws is thus Nsv x T x Nmixtures
% zdraws      = bsxfun(@minus, logy2 - h - outlierlog2Draws, KSCt.mean) ./ KSCt.vol;
zdraws      = (logy2 - h - SVtScalelog2Draws - KSCt.mean) ./ KSCt.vol;

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
obs   = logy2 - KSC.mean(kai2States) - SVtScalelog2Draws;
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


%% construct SV
SV = exp((h + SVtScalelog2Draws) / 2);

