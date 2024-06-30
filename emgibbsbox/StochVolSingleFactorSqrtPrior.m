function [h, lambda, lambdaResid, hbar, kai2States] = StochVolSingleFactorSqrtPrior(logy2, h, beta, lambdavar, Eh0, sqrtVh0, KSC, KSCt, Nsv, T, rndStream)
% StochVolSingleFactor ...
%

%   Coded by  Elmar Mertens, em@elmarmertens.com


%% VERSION INFO
% AUTHOR    : Elmar Mertens
% $DATE     : 28-Aug-2009 12:07:01 $
% $Revision : 1.00 $
% DEVELOPED : 7.7.0.471 (R2008b)
% FILENAME  : StochVolKSC.m.m

%% draw mixture states
% zdraws are standardized draws for each component of the normal mixture 
% zdraws is thus Nsv x T x 7 
zdraws      = (logy2 - h - KSCt.mean) ./ KSCt.vol;

% construct CDF
% factor of sqrt(2 * pi) can be ommitted for kernel
pdfKernel           = KSCt.pdf ./ KSCt.vol .* exp(-.5 * zdraws.^2); 
cdf                 = cumsum(pdfKernel, 3);                % integrate
cdf(:,:,1:end-1)    = cdf(:,:,1:end-1) ./ cdf(:,:,end); 
cdf(:,:,end)        = 1;    % normalize

% draw states
kai2States  = sum(rand(rndStream, Nsv, T) > cdf, 3) + 1;
obs         = logy2 - KSC.mean(kai2States);

%% Single-Factor State Space
A      = eye(Nsv+1);
B      = zeros(Nsv+1,1);
B(1,1) = sqrt(lambdavar);
C      = cat(2, beta, eye(Nsv));
sqrtR  = zeros(Nsv,Nsv,T);
for n = 1 : Nsv
    sqrtR(n,n,:) = KSC.vol(kai2States(n,:));
end

EX0     = cat(1, 0, Eh0);
sqrtVX0 = blkdiag(0, sqrtVh0);

[X, Xshock, X0] = a2b2c2DisturbanceSmoothingSampler1draw(A, B, C, obs, EX0, sqrtVX0, ...
    sqrtR, rndStream); 

lambda      = X(1,:);
lambdaResid = Xshock(1,:);
hbar        = X0(2:end);
h           = beta * lambda + hbar;



