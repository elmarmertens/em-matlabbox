function [Xdraw, Xshock, X0draw, Xhat, P] = commonAR1noisePrecisionBasedSampler(Y, Ny, T, rhoSTATE, volSTATE, volNOISE, Ndraws, rndStream)
% precisionBasedSampler for common AR1 state (plus noise)
% computes smoothed kalman states using the stacked approach of Chan and Jeliazkov
% Chan and Jeliazkov
%
%   ...

% assumes rhoSTATE and volSTATE are scalars, and volNOISE is (Ny x T) x 1 vector (i.e. no correlation within X and Y)
% X0 is zero but sqrtV0 is non-zero (allowing stochastic initial condition)
% Note: mean of AR1 is zero (initial condition is free)

if nargin < 7 || isempty(Ndraws)
    Ndraws = 1;
end
if nargin < 8
    rndStream = getDefaultStream;
end

%% read parameters
Y        = Y(:);
NyT      = Ny * T;
Nx       = 1;

%% construct stacked system
XX0 = sparse(T+1, 1);

% AA
rowndx = [1 : T+1, Nx + (1 : T)];
colndx = [1 : T+1, 1 : T];
values = [ones(1,T+1), -rhoSTATE * ones(1, T)];
AA     = sparse(rowndx, colndx, values);

% CC
rowndx = 1 : NyT;
colndx = repmat(1:T, Ny, 1);
colndx = colndx(:) + 1;
values = ones(1, NyT);
CC     = sparse(rowndx, colndx, values);


if isscalar(volSTATE)
    sqrtSIGMA    = volSTATE * speye(T+1);
    sqrtSIGMA(1) = 10 * volSTATE;
else
    volSTATE    = [10 * volSTATE(1); volSTATE(:)]; % inflate initial state variance
    sqrtSIGMA   = sparse(1:T+1, 1:T+1, volSTATE(:));
end
sqrtOMEGA    = sparse(1:NyT, 1:NyT, volNOISE(:)');

%% set up  stacked system


AAtilde            = sqrtSIGMA \ AA;
% XX0tilde           = sqrtSIGMA \ XX0; % note: could exploit XX0=XX0tilde=0

CCtilde            = sqrtOMEGA \ CC;
Ytilde             = sqrtOMEGA \ Y;

P                   = AAtilde' * AAtilde + (CCtilde' * CCtilde);
[sqrtP, flag]       = chol(P, 'lower');

if flag > 0
    error('P not posdf, using QR instead')
    % via qr -- much slower
    M = [AAtilde; CCtilde]; %#ok<UNRCH>
    m = size(M,2);
    [~, R] = qr(M);
    sqrtP = R(1:m,1:m)';
    % checkdiff(sqrtP * sqrtP', sqrtP2 * sqrtP2');
end

% sqrtPXhat   = sqrtP \ (AAtilde' * XX0tilde + CCtilde' * Ytilde); % note: could exploit XX0=XX0tilde=0
sqrtPXhat   = sqrtP \ (CCtilde' * Ytilde); % note: exploiting XX0=XX0tilde=0

Zdraw        = randn(rndStream, T+1, Ndraws);
XXdraw       = (sqrtP') \ (sqrtPXhat + Zdraw);
Xdraw        = XXdraw(2:end);

if nargout > 1
    Xshock = AA * XXdraw;
    Xshock = reshape(Xshock(2:end), Nx, T, Ndraws); % ignore initial value
end
Xdraw        = reshape(Xdraw, Nx, T, Ndraws);

if nargout > 2
    X0draw       = XXdraw(1);
end

if nargout > 3
    Xhat        = (sqrtP') \ sqrtPXhat;
    Xhat        = reshape(Xhat(2:end), Nx, T);
end

