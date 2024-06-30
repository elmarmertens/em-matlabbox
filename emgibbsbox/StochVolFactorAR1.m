function [h, lambda, lambdaShock, hbar, kai2States] = StochVolFactorAR1(logy2, h, beta, lambdarho, lambdavar, Eh0, Vh0, KSC, KSCt, Nsv, T, rndStream)
% StochVolSingleFactor ...
%

%   Coded by  Elmar Mertens, em@elmarmertens.com


% todo: currently ignoring lambdarho and using RW instead of AR1

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

%% Single-Factor State Space, with AR(1) SV factor
Nfactors = size(beta,2);

obs                      = logy2 - KSC.mean(kai2States);
A                        = eye(Nfactors+Nsv);
A(1:Nfactors,1:Nfactors) = diag(lambdarho);

B      = cat(1, chol(lambdavar, 'lower'), zeros(Nsv,Nfactors));
C      = cat(2, beta, eye(Nsv));
sqrtR  = zeros(Nsv,Nsv,T);
for n = 1 : Nsv
    sqrtR(n,n,:) = KSC.vol(kai2States(n,:));
end

EX0     = cat(1, zeros(Nfactors,1), Eh0);
sqrtVX0 = diag(cat(1, zeros(Nfactors,1), sqrt(Vh0)));

[X, Xshock, X0] = a2b2c2DisturbanceSmoothingSampler1draw(A, B, C, obs, EX0, sqrtVX0, ...
    sqrtR, rndStream); 

lambda      = X(1:Nfactors,:);
lambdaShock = Xshock(1:Nfactors,:);
hbar        = X0(Nfactors+1:end);
h           = beta * lambda + hbar;



