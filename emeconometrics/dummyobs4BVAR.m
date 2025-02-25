function [Ystar, Xstar] = dummyobs4BVAR(Ny, p, Nx, ndxExo, lambda, priorMean, dof0, SIGMA0)
% BAYESREGDUMMYOBS constructs dummy obs for Bayesian Regression with Conjugate Prior
%
% usage [Ystar, Xstar] = dummyobs4BVAR(Ny, p, Nx, ndxExo, lambda, priorMean, dof0, SIGMA0)
%
% Note: ndxExo can be a constant intercept or other deterministic dummies
%
% defaults ndxExo=[], priorMean=0, dof=0 (SIGMA0 irrelevant)
%
% see also BVARjeffries, BVARdraws

%% VERSION INFO
% AUTHOR    : Elmar Mertens
% $DATE     : 04-Nov-2021 15:40:03 $
% $Revision : 1.00 $
% DEVELOPED : 9.11.0.1769968 (R2021b)
% FILENAME  : dummyobs4BVAR.m

%% process arguments

if nargin < 4
    ndxExo = [];
end

if nargin < 5 || isempty(lambda)
    % hyperparameters (in Sims notation)
    lambda0     = 1/10;
    lambda1     = 5; % overall shrinkage, corresponds to 1/lambda1^2 = .2^2 in CCMM or RATS
    lambda2     = 1; % decay, corresponds to lambda2*2=2 in CCMM or RATS
    % lambda = [lambda0 lambda1 lambda2];
else
    lambda0 = lambda(1);
    lambda1 = lambda(2);
    lambda2 = lambda(3);
end

if nargin < 6 || isempty(priorMean)
    priorMean = 0;
end

if nargin < 7 || isempty(dof0)
    dof0 = 0;
end

if nargin < 8 || isempty(SIGMA0)
    SIGMA0 = eye(Ny);
end

if isscalar(priorMean)
    priorMean = repmat(priorMean, 1, Ny);
end

if islogical(ndxExo)
    ndxExo = find(ndxExo);
end

if isempty(ndxExo)
    ndxLags = 1 : Nx;
else
    ndxLags = ~ismember(1 : Nx, ndxExo);
end

%% construct dummy obs for Bayesian prior

% Set up dummy variables
Xstar = [];
Ystar = [];

sqrtSIGMA0 = chol(SIGMA0)';

if dof0 > 0
    if dof0 < Ny
        error('IW prior must have at least Ny dof')
    end
    thisYstar         = zeros(dof0, Ny);
    thisYstar(1:Ny,:) = sqrtSIGMA0 * dof0; 
    % note: 
    % - the above is numerically equivalent to stacking sqrtSIGMA0 dof0 / Ny times, 
    % - but can also handle cases where dof0 is not a multiple of Ny
    % - with dof0>=Ny there are potentially rows padded with zeros at bottom of thisYstar (needed to compute correct dof = Tstar-k)
    
    
    Ystar             = cat(1, Ystar, thisYstar);
    Xstar             = cat(1, Xstar, zeros(dof0, Nx));
end


% persistence priors
thisYstar = zeros(Ny * p, Ny);
thisXstar = zeros(Ny * p, Nx);

% build thisXstar via kronecker
LAMBDA               = lambda1 .* ((1:p) .^ lambda2);
thisXstar(:,ndxLags) = kron(diag(LAMBDA), sqrtSIGMA0);


% unit root on own persistence
thisYstar(1:Ny,:) = lambda1 * sqrtSIGMA0 .* diag(priorMean);

Ystar = cat(1, Ystar, thisYstar);
Xstar = cat(1, Xstar, thisXstar);

% intercept prior
if lambda0 > 0
    Nexo               = length(ndxExo);
    constXstar         = zeros(Nexo,Nx);
    for n = 1 : Nexo
        constXstar(n,ndxExo(n)) = lambda0;
    end

    Ystar = cat(1, Ystar, zeros(Nexo,Ny));
    Xstar = cat(1, Xstar, constXstar);
end