function Xdraws = stateABnanDraw1(A, B, ndxY, Ydata, yNaNndx, X0, sqrtSigma, rndStream)
% STATEABCDRAW
% ....
% supposes X0 is deterministically given

% note: rows of C with yNaNndx == true must be zero (no checks!)

%   Coded by  Elmar Mertens, em@elmarmertens.com


%% parse inputs
Nx                = size(A, 1);
[Ny, T]           = size(Ydata);
Nw                = size(B,2);

if nargin < 8 || isempty(rndStream)
    rndStream = getDefaultStream;
end



%% init Variables and allocate memory
I          = eye(Nx);

yDataNdx          = ~yNaNndx;
obsndx = false(Nx,T);
for t = 1 : T
    obsndx(ndxY(yDataNdx(:,t)),t) = true;
end

%% allocate memory
[Sigmattm1, ImKC]           = deal(zeros(Nx, Nx, T));
invSigmaYttm1               = zeros(Ny, Ny, T);
Ytilde                      = zeros(Ny, T);
[XtT, Xttm1, Xplus]         = deal(zeros(Nx, T));


%% generate plus data
wplus   = randn(rndStream, Nw, T);

%% Forward Loop: Kalman Forecasts
Sigmatt = zeros(Nx,Nx);
Xtt     = zeros(Nx,1); % mean of differences between actual and plus (and thus zeros)

disturbanceplus  = zeros(Nx, T);

for t = 1 : T
    
    % "plus" States and priors
    Bsv                     = B * diag(sqrtSigma(:,t));
    disturbanceplus(:,t)    = Bsv * wplus(:,t);
    BSigmaB                 = Bsv * Bsv';
    
    if t == 1
        Xlagplus = X0;
    else
        Xlagplus = Xplus(:,t-1);
    end
    
    Xplus(:,t)              = A * Xlagplus + disturbanceplus(:,t);
    
    % priors
    Sigmattm1(:,:,t)        = A * Sigmatt * A' + BSigmaB;
    Xttm1(:,t)              = A * Xtt;
    
    % observed innovation
    Yplus            = Xplus(obsndx(:,t),t);
    SigmaYttm1       = Sigmattm1(obsndx(:,t),obsndx(:,t),t);
    
    Ytilde(yDataNdx(:,t),t)      = Ydata(yDataNdx(:,t),t) - Yplus - Xttm1(obsndx(:,t),t);
    
    % Block Inverse of Y-VCV, accounting for missing obs (with zero VCV)
    invSigmaYttm1(yDataNdx(:,t),yDataNdx(:,t),t) = eye(sum(yDataNdx(:,t))) / SigmaYttm1;
    
     
    % Kalman Gain
    K                       = Sigmattm1(:,obsndx(:,t),t) * invSigmaYttm1(yDataNdx(:,t),yDataNdx(:,t),t);

    % posteriors
    Xtt                     = Xttm1(:,t) + K * Ytilde(yDataNdx(:,t),t);

    KC                      = zeros(Nx);
    KC(:,obsndx(:,t))       = K;
    ImKC(:,:,t)             = I-KC;
    Sigmatt                 = ImKC(:,:,t) * Sigmattm1(:,:,t) * ImKC(:,:,t)'; % Joseph form for better numerical stability
    
end

%% Backward Loop: Disturbance Smoother
XtT(:,T)   = Xtt;

t = T;
Yscaled          = invSigmaYttm1(yDataNdx(:,t),yDataNdx(:,t),t) * Ytilde(yDataNdx(:,t),t);
StT              = zeros(Nx,1);
StT(obsndx(:,t)) = Yscaled;

for t = (T-1) : -1 : 1
    
    Atilde          = A * ImKC(:,:,t);
    
    Yscaled             = invSigmaYttm1(yDataNdx(:,t),yDataNdx(:,t),t) * Ytilde(yDataNdx(:,t),t);
    stilde              = zeros(Nx,1);
    stilde(obsndx(:,t)) = Yscaled;
    StT                 = Atilde' * StT + stilde;
   
    XtT(:,t)        = Xttm1(:,t) + Sigmattm1(:,:,t) * StT;

end

%% sample everything together (and reorder output dimensions)
Xdraws  = XtT + Xplus;



