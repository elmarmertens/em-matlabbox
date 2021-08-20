function Xdraws = stateABnanDraws(A, B, ndxY, Ydata, yNaNndx, X0, sqrtSigma, Ndraws, rndStream)
% STATEABCDRAW
% ....
% supposes X0 is deterministically given

% note: rows of C with yNaNndx == true must be zero (no checks!)

%   Coded by  Elmar Mertens, em@elmarmertens.com


%% parse inputs
Nx                = size(A, 1);
[Ny, T]           = size(Ydata);
Ydata             = permute(Ydata, [1 3 2]); 
Nw                = size(B,2);

if nargin < 8 || isempty(Ndraws)
    Ndraws = 1;
end
if nargin < 9 || isempty(rndStream)
    rndStream = getDefaultStream;
end



%% init Variables and allocate memory
I                 = eye(Nx);

yDataNdx          = ~yNaNndx;
obsndx = false(Nx,T);
for t = 1 : T
    obsndx(ndxY(yDataNdx(:,t)),t) = true;
end

CC         = zeros(Ny,Nx);
CC(:,ndxY) = eye(Ny);


%% allocate memory
[Sigmattm1, ImKC]           = deal(zeros(Nx, Nx, T));
invSigmaYttm1               = zeros(Ny, Ny, T);
Ytilde                      = zeros(Ny, 1, T);
Yplustilde                  = zeros(Ny, Ndraws, T);
[XtT, Xttm1]                = deal(zeros(Nx, 1, T));
[Xplus, XplustT, Xplusttm1] = deal(zeros(Nx, Ndraws, T));

Sigmatt          = zeros(Nx,Nx);
Xtt              = X0;
Xplustt          = repmat(X0, 1, Ndraws);
disturbanceplus  = zeros(Nx, Ndraws, T);


%% draw random numbers for plus data

wplus   = randn(rndStream, Nw, Ndraws, T);

%% Forward Loop: Kalman Forecasts

for t = 1 : T
    
    % "plus" States and priors
    Bsv                     = B * diag(sqrtSigma(:,t));
    disturbanceplus(:,:,t)  = Bsv * wplus(:,:,t);
    BSigmaB                 = Bsv * Bsv';
    
    if t == 1
        Xlagplus = X0;
    else
        Xlagplus = Xplus(:,:,t-1);
    end
    Xplus(:,:,t)           = A * Xlagplus + disturbanceplus(:,:,t);

    % priors
    Sigmattm1(:,:,t)        = A * Sigmatt * A' + BSigmaB;
    Xttm1(:,:,t)            = A * Xtt;
    Xplusttm1(:,:,t)        = A * Xplustt;
    
    
    % observed innovation
    Yplus            = Xplus(obsndx(:,t),:,t);
    SigmaYttm1       = Sigmattm1(obsndx(:,t),obsndx(:,t),t);
    
    %     that = C(:,:,t) * Xplus(:,t);
    %     checkdiff(Yplus, that(yDataNdx(:,t)));
    %     that =  C(:,:,t) * Sigmattm1(:,:,t) * C(:,:,t)';
    %     checkdiff(SigmaYttm1, that(yDataNdx(:,t),yDataNdx(:,t)));
    
    Ytilde(yDataNdx(:,t),:,t)      = Ydata(yDataNdx(:,t),:,t) - Xttm1(obsndx(:,t),:,t);
    Yplustilde(yDataNdx(:,t),:,t)  = Yplus - Xplusttm1(obsndx(:,t),:,t);
    
    % Block Inverse of Y-VCV, accounting for missing obs (with zero VCV)
    invSigmaYttm1(yDataNdx(:,t),yDataNdx(:,t),t) = eye(sum(yDataNdx(:,t))) / SigmaYttm1;
    
     
    % Kalman Gain
    K                       = Sigmattm1(:,obsndx(:,t),t) * invSigmaYttm1(yDataNdx(:,t),yDataNdx(:,t),t);
    ImKC(:,:,t)             = I - K * CC(yDataNdx(:,t),:);    
    
    % posteriors
    Sigmatt                 = ImKC(:,:,t) * Sigmattm1(:,:,t) * ImKC(:,:,t)'; % Joseph form for better numerical stability
    Xtt                     = Xttm1(:,:,t)     + K * Ytilde(yDataNdx(:,t),:,t);
    Xplustt                 = Xplusttm1(:,:,t) + K * Yplustilde(yDataNdx(:,t),:,t);
    
end

%% Backward Loop: Disturbance Smoother

t = T;
XplustT(:,:,t)  = Xplustt;
XtT(:,:,t)      = Xtt;
thisC           = CC(yDataNdx(:,t),:);
StT             = thisC' * (invSigmaYttm1(yDataNdx(:,t),yDataNdx(:,t),t) * Ytilde(yDataNdx(:,t),:,t));
SplustT         = thisC' * (invSigmaYttm1(yDataNdx(:,t),yDataNdx(:,t),t) * Yplustilde(yDataNdx(:,t),:,t));


for t = (T-1) : -1 : 1
    
    Atilde          = A * ImKC(:,:,t);
    thisC           = CC(yDataNdx(:,t),:);
    
    StT             = Atilde' * StT + thisC' * (invSigmaYttm1(yDataNdx(:,t),yDataNdx(:,t),t) * Ytilde(yDataNdx(:,t),:,t));
    XtT(:,:,t)      = Xttm1(:,:,t) + Sigmattm1(:,:,t) * StT;
    
    SplustT         = Atilde' * SplustT + thisC' * (invSigmaYttm1(yDataNdx(:,t),yDataNdx(:,t),t) * Yplustilde(yDataNdx(:,t),:,t));
    XplustT(:,:,t)  = Xplusttm1(:,:,t) + Sigmattm1(:,:,t) * SplustT;
end

%% sample everything together (and reorder output dimensions)
Xdraws  = bsxfun(@minus, XtT, XplustT) + Xplus;
Xdraws  = permute(Xdraws, [1 3 2]);



