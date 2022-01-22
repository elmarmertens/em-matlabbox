function [Xdraw, Xshock, Xhat, P] = commonAR1noisePrecisionBasedSampler(Y, Ny, T, rhoSTATE, volSTATE, volNOISE, Ndraws, rndStream)
% abcdPrecisionBasedSampler computes smoothed kalman states using the stacked approach of 
% Chan and Jeliazkov
%  
%   ... 

% assumes rhoSTATE and volSTATE are scalars, and volNOISE is (Ny x T) x 1 vector (i.e. no correlation within X and Y)
% X0 and sqrtV0 both zero, so that X0=0 is deterministic

if nargin < 7 || isempty(Ndraws)
    Ndraws = 1;
end
if nargin < 8
    rndStream = getDefaultStream;
end

%% read parameters
Y        = Y(:);
NyT      = Ny * T;
Nx       = 1; % todo: hardcode
NxT      = Nx * T; % obsolete but for better readability
NxTm1    = Nx * (T - 1); % obsolete but for better readability

%% construct stacked system
XX0 = sparse(T, 1);

% AA
rowndx = [1 : NxT, Nx + 1 : NxT];
colndx = [1 : NxT, 1 : NxTm1];
values = [ones(1,NxT), -rhoSTATE * ones(1, T-1)];
AA     = sparse(rowndx, colndx, values);

% CC
rowndx = 1 : NyT;
colndx = repmat(1:T, Ny, 1);
colndx = colndx(:);
values = ones(1, NyT);
CC     = sparse(rowndx, colndx, values);



sqrtSIGMA   = volSTATE * speye(NxT); % Note: X(t=1) consists solely of innovation since X(t=0)=0
sqrtOMEGA   = sparse(1:NyT, 1:NyT, volNOISE(:)');

%% set up  stacked system


AAtilde            = sqrtSIGMA \ AA;
XX0tilde           = XX0; % sqrtSIGMA \ XX0; since XX0=0

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

sqrtPXhat   = sqrtP \ (AAtilde' * XX0tilde + CCtilde' * Ytilde); 

Zdraw        = randn(rndStream, Nx * T, Ndraws);
Xdraw        = (sqrtP') \ (sqrtPXhat + Zdraw);

if nargout > 1
    Xshock = AA * Xdraw;
    Xshock = reshape(Xshock, Nx, T, Ndraws);
end

Xdraw        = reshape(Xdraw, Nx, T, Ndraws);

if nargout > 2
    Xhat        = (sqrtP') \ sqrtPXhat;
    Xhat        = reshape(Xhat, Nx, T);
end

