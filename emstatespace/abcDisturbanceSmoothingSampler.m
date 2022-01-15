function [Xdraws, disturbanceDraws, X0draws, noiseDraws] = ...
    abcDisturbanceSmoothingSampler(A, B, C, Ydata, X00, cholSigma00, Ndraws, ...
    sqrtR, sqrtSigma, rndStream)
% ABCDISTURBANCESMOOTHINGSAMPLER: samples Ndraws from ABC-R system (w/optional SV)
% 
% usage: [Xdraws, disturbanceDraws, X0draws, noiseDraws] = ...
%    abcDisturbanceSmoothingSampler(A, B, C, Ydata, X00, cholSigma00, Ndraws, ...
%    sqrtR, sqrtSigma, rndStream)
% where: sqrtSigma and rndStream are optional
%    
% see also abcDisturbanceSmoothingSampler1draw

%   Coded by  Elmar Mertens, em@elmarmertens.com



%% parse inputs
Nx                = size(A, 1);
[Ny, T]           = size(Ydata);
Ydata             = permute(Ydata, [1 3 2]);

Nw                = size(B,2);

if nargin < 7 || isempty(Ndraws)
    Ndraws = 1;
end

if nargin < 8
    sqrtR = [];
end

if nargin < 9
    sqrtSigma = [];
end

if nargin < 10
    rndStream = RandStream.getGlobalStream;
end

%% init Variables and allocate memory
if ~isempty(sqrtR) && ismatrix(sqrtR)
    sqrtR = repmat(sqrtR, [1 1 T]);
end
if ~isempty(sqrtSigma) && isvector(sqrtSigma)
    sqrtSigma = repmat(sqrtSigma(:), [1 T]);
end

if ismatrix(A)
    A = repmat(A, [1 1 T]);
end
if ismatrix(B)
    B = repmat(B, [1 1 T]);
end
if ismatrix(C)
    C = repmat(C, [1 1 T]);
end

I                 = eye(Nx);
Iy                = eye(Ny);

%% allocate memory
[Sigmattm1, ImKC]           = deal(zeros(Nx, Nx, T));
invSigmaYttm1               = zeros(Ny, Ny, T);
Ytilde                      = zeros(Ny, 1, T);
Yplustilde                  = zeros(Ny, Ndraws, T);
[XtT, Xttm1]                = deal(zeros(Nx, 1, T));
[Xplus, XplustT, Xplusttm1] = deal(zeros(Nx, Ndraws, T));


%% generate plus data

wplus   = randn(rndStream, Nw, Ndraws, T);
if ~isempty(sqrtR)
    Nmeasurmentnoise = size(sqrtR, 2);
    eplus    = randn(rndStream, Nmeasurmentnoise, Ndraws, T);
end
X0plus = X00 + cholSigma00 * randn(rndStream, Nx, Ndraws);

%% Forward Loop: Kalman Forecasts
[Sigma00, Sigmatt] = deal(cholSigma00 * cholSigma00');
Xtt     = X00;
Xplustt = repmat(X00, 1, Ndraws);
BSigmaB = zeros(Nx, Nx, T);

disturbanceplus  = zeros(Nx, Ndraws, T);
if isempty(sqrtR)
    noiseplus        = [];
else
    noiseplus        = zeros(Ny, Ndraws, T);
end

for t = 1 : T

    % "plus" States and priors
    if isempty(sqrtSigma)
        disturbanceplus(:,:,t)  = B(:,:,t) * wplus(:,:,t);
        BSigmaB(:,:,t)          = B(:,:,t) * B(:,:,t)';
    else
        disturbanceplus(:,:,t)  = B(:,:,t) * diag(sqrtSigma(:,t)) * wplus(:,:,t);
        BSigmaB(:,:,t)          = B(:,:,t) * diag(sqrtSigma(:,t).^2) * B(:,:,t)';
    end

    if t == 1
        Xplus(:,:,t)      = A(:,:,t) * X0plus + disturbanceplus(:,:,t);
    else
        Xplus(:,:,t)      = A(:,:,t) * Xplus(:,:,t-1) + disturbanceplus(:,:,t);
    end

    % priors
    Sigmattm1(:,:,t)     = A(:,:,t) * Sigmatt * A(:,:,t)' + BSigmaB(:,:,t);
    Xttm1(:,:,t)         = A(:,:,t) * Xtt;
    Xplusttm1(:,:,t)     = A(:,:,t) * Xplustt;


    % observed innovation
    if isempty(sqrtR)

        Yplus            = C(:,:,t) * Xplus(:,:,t);
        SigmaYttm1       = C(:,:,t) * Sigmattm1(:,:,t) * C(:,:,t)';

    else

        noiseplus(:,:,t)  = sqrtR(:,:,t) * eplus(:,:,t); %#ok<AGROW>
        Yplus             = C(:,:,t) * Xplus(:,:,t) + noiseplus(:,:,t);
        SigmaYttm1        = C(:,:,t) * Sigmattm1(:,:,t) * C(:,:,t)' + sqrtR(:,:,t) * sqrtR(:,:,t)';

    end

    Ytilde(:,:,t)         = Ydata(:,:,t)    - C(:,:,t) * Xttm1(:,:,t);
    Yplustilde(:,:,t)     = Yplus           - C(:,:,t) * Xplusttm1(:,:,t);

    % Block Inverse of Y-VCV, accounting for missing obs (with zero VCV)
    invSigmaYttm1(:,:,t) = Iy / SigmaYttm1;

    % Kalman Gain
    K                       = (Sigmattm1(:,:,t) * C(:,:,t)') * invSigmaYttm1(:,:,t);
    ImKC(:,:,t)             = I - K * C(:,:,t);

    % posteriors
    Sigmatt                 = ImKC(:,:,t) * Sigmattm1(:,:,t);

    Xtt                     = Xttm1(:,:,t) + K * Ytilde(:,:,t);
    Xplustt                 = Xplusttm1(:,:,t) + K * Yplustilde(:,:,t);

end

%% Backward Loop: Disturbance Smoother
XplustT(:,:,T)  = Xplustt;
XtT(:,T)        = Xtt;

StT             = C(:,:,T)' * (invSigmaYttm1(:,:,T) * Ytilde(:,T));
SplustT         = C(:,:,T)' * (invSigmaYttm1(:,:,T) * Yplustilde(:,:,T));

if nargout > 1
    disturbancetT               = zeros(Nx, 1, T);
    disturbancetT(:,T)          = BSigmaB(:,:,T) * StT;
    disturbanceplustT           = zeros(Nx, Ndraws, T);
    disturbanceplustT(:,:,T)    = BSigmaB(:,:,T) * SplustT;
else
    disturbancetT        = [];
    disturbanceplustT    = [];
end

if nargout > 2 && ~isempty(sqrtR)
    noisetT              = zeros(Ny, 1, T);
    noisetT(:,:,T)       = Ytilde(:,:,T) - C(:,:,T) * (XtT(:,:,T) - Xttm1(:,:,T));
    noiseplustT          = zeros(Ny, Ndraws, T);
    noiseplustT(:,:,T)   = Yplustilde(:,:,T) - C(:,:,T) * (XplustT(:,:,T) - Xplusttm1(:,:,T));
else
    noisetT      = [];
    noiseplustT  = [];
end


for t = (T-1) : -1 : 1
    Atilde      = A(:,:,t+1) * ImKC(:,:,t);

    StT         = Atilde' * StT + ...
        C(:,:,t)' * (invSigmaYttm1(:,:,t) * Ytilde(:,:,t));
    XtT(:,:,t)    = Xttm1(:,:,t) + Sigmattm1(:,:,t) * StT;

    SplustT     = Atilde' * SplustT + ...
        C(:,:,t)' * (invSigmaYttm1(:,:,t) * Yplustilde(:,:,t));
    XplustT(:,:,t)  = Xplusttm1(:,:,t) + Sigmattm1(:,:,t) * SplustT;


    if ~isempty(disturbancetT)
        disturbancetT(:,:,t)        = BSigmaB(:,:,t) * StT;
        disturbanceplustT(:,:,t)    = BSigmaB(:,:,t) * SplustT;
    end
    if ~isempty(noisetT)
        noisetT(:,:,t)       = Ytilde(:,:,t) - C(:,:,t) * (XtT(:,:,t) - Xttm1(:,:,t));
        noiseplustT(:,:,t)   = Yplustilde(:,:,t) - C(:,:,t) * (XplustT(:,:,t) - Xplusttm1(:,:,t));
    end

end

%% sample everything together (and reorder output dimensions)
Xdraws  = permute(bsxfun(@minus, XtT, XplustT) + Xplus, [1 3 2]);

if nargout > 1

    disturbanceDraws  = permute(bsxfun(@minus, disturbancetT, disturbanceplustT) ...
        + disturbanceplus, [1 3 2]);
    if any(isnan(disturbanceDraws(:)))
        error('em:debug','wDraws contain NaNs ...');
    end

    if nargout > 2

        X0T      = X00 + Sigma00 * A(:,:,1)' * StT;
        Xplus0T  = bsxfun(@plus, X00, Sigma00 * A(:,:,1)' * SplustT);

        X0draws  = bsxfun(@minus, X0T, Xplus0T) + X0plus;


        if nargout > 3
            noiseDraws  = permute(bsxfun(@minus, noisetT, noiseplustT) ...
                + noiseplus, [1 3 2]);
        end
    end
end
