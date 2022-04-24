function [XtT, SigmatT, Sigmatt, Sigmattm1] = abcDisturbanceSmootherNaN(A, B, C, Ydata, yNaNndx, X0, Sigma0)
% See also abcStateSmoother

%   Coded by  Elmar Mertens, em@elmarmertens.com

%% init Variables and allocate memory
Nx      = size(A, 1);
[Ny, T] = size(Ydata);
% Nw      = size(B, 2);
I       = eye(Nx);


[Sigmattm1, Sigmatt, KC]   = deal(zeros(Nx, Nx, T));
K                          = zeros(Nx, Ny, T);
invSigmaYttm1              = zeros(Ny, Ny, T);
Ytilde                     = zeros(Ny, T);
[XtT, Xtt, Xttm1, StT]     = deal(zeros(Nx, T));

yDataNdx = ~yNaNndx;


%% Forward Loop: Kalman Forecasts
for t = 1 : T
    
    % setup priors
    if t == 1
        Sigmattm1(:,:,1)    = Sigma0;
        Xttm1(:,1)          = X0;
    else
        Sigmattm1(:,:,t)    = A(:,:,t) *  Sigmatt(:,:,t-1) * A(:,:,t)' + B(:,:,t) * B(:,:,t)';
        Xttm1(:,t)          = A(:,:,t) * Xtt(:,t-1);
    end
    
    % observed innovation
    Ytilde(:,t)         = Ydata(:,t) - C(:,:,t) * Xttm1(:,t);
    % Kalman Gain
    SigmaYttm1          = C(:,:,t) * Sigmattm1(:,:,t) * C(:,:,t)';
    % Block Inverse of Y-VCV, accounting for missing obs (with zero VCV)
    invSigmaYttm1(yDataNdx(:,t),yDataNdx(:,t),t)  = eye(sum(yDataNdx(:,t))) ...
        / SigmaYttm1(yDataNdx(:,t),yDataNdx(:,t));
    K(:,:,t)            = (Sigmattm1(:,:,t) * C(:,:,t)') * invSigmaYttm1(:,:,t);
    KC(:,:,t)           = K(:,:,t) * C(:,:,t);
    % posteriors
    Sigmatt(:,:,t)      = (I - KC(:,:,t)) * Sigmattm1(:,:,t);
    Xtt(:,t)            = Xttm1(:,t) + K(:,:,t) * Ytilde(:,t);

    %     checkdiff(Sigmatt(:,:,t), ...
    %         (I - K(:,:,t) * C(:,:,t)) * Sigmattm1(:,:,t) * (I - K(:,:,t) * C(:,:,t))');
    %     checkdiff(Sigmatt(:,:,t), ...
    %         Sigmattm1(:,:,t) * (I - K(:,:,t) * C(:,:,t))');
 
end


%% Backward Loop: Disturbance Smoother
StT(:,T) = C(:,:,T)' * (invSigmaYttm1(:,:,T) *  Ytilde(:,T)); 
XtT(:,T) = Xtt(:,T);
if nargout > 1
   SigmatT        = zeros(Nx, Nx, T);
   N              = C(:,:,T)' * invSigmaYttm1(:,:,T) * C(:,:,T);
   SigmatT(:,:,T) = Sigmatt(:,:,T);
   checkdiff(SigmatT(:,:,T), Sigmattm1(:,:,T) - Sigmattm1(:,:,T)  * N * Sigmattm1(:,:,T));
end

for t = T -1 : -1 : 1
    Atilde   = A(:,:,t+1) * (I - KC(:,:,t));
    StT(:,t) = C(:,:,t)' * (invSigmaYttm1(:,:,t) *  Ytilde(:,t)) + Atilde' * StT(:,t+1);
    XtT(:,t) = Xttm1(:,t) + Sigmattm1(:,:,t) * StT(:,t);
    if nargout > 1
        N               = C(:,:,t)' * (invSigmaYttm1(:,:,t) *  C(:,:,t)) + Atilde' * N * Atilde;
        SigmatT(:,:,t)  = Sigmattm1(:,:,t) - Sigmattm1(:,:,t)  * N * Sigmattm1(:,:,t);
    end
end
