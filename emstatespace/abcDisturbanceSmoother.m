function [XtT, SigmatT, Xtt, Sigmatt, Sigmattm1] = abcDisturbanceSmoother(A, B, C, Ydata, X0, Sigma0)
% See also abcStateSmootherNAN

%   Coded by  Elmar Mertens, em@elmarmertens.com

%% init Variables and allocate memory
Nx      = size(A, 1);
[Ny, T] = size(Ydata);
I       = eye(Nx);

if ismatrix(A)
    A = repmat(A, [1 1 T]);
end
if ismatrix(B)
    B = repmat(B, [1 1 T]);
end
if ismatrix(C)
    C = repmat(C, [1 1 T]);
end

[Sigmattm1, Sigmatt, KC]   = deal(zeros(Nx, Nx, T));
K           = zeros(Nx, Ny, T);
SigmaYttm1  = zeros(Ny, Ny, T);
Ytilde      = zeros(Ny, T);
[XtT, Xtt, Xttm1, StT] = deal(zeros(Nx, T));

%% Forward Loop: Kalman Forecasts
for t = 1 : T
    
    % setup priors
    if t == 1
        Sigmattm1(:,:,1)    = A(:,:,t) *  Sigma0 * A(:,:,t)' + B(:,:,t) * B(:,:,t)';
        Xttm1(:,1)          = A(:,:,t) * X0;
    else
        Sigmattm1(:,:,t)    = A(:,:,t) *  Sigmatt(:,:,t-1) * A(:,:,t)' + B(:,:,t) * B(:,:,t)';
        Xttm1(:,t)          = A(:,:,t) * Xtt(:,t-1);
    end
    
    % observed innovation
    Ytilde(:,t)         = Ydata(:,t) - C(:,:,t) * Xttm1(:,t);
    % Kalman Gain
    SigmaYttm1(:,:,t)   = C(:,:,t) * Sigmattm1(:,:,t) * C(:,:,t)';
    K(:,:,t)            = (Sigmattm1(:,:,t) * C(:,:,t)') / SigmaYttm1(:,:,t);
    KC(:,:,t)           = K(:,:,t) * C(:,:,t);
    % posteriors
    Sigmatt(:,:,t)      = (I - KC(:,:,t)) * Sigmattm1(:,:,t) * (I - K(:,:,t) * C(:,:,t))'; % Joseph Form for numerical stability
    Xtt(:,t)            = Xttm1(:,t) + K(:,:,t) * Ytilde(:,t);

    %     checkdiff(Sigmatt(:,:,t), ...
    %         (I - K(:,:,t) * C(:,:,t)) * Sigmattm1(:,:,t) * (I - K(:,:,t) * C(:,:,t))');
    %     checkdiff(Sigmatt(:,:,t), ...
    %         Sigmattm1(:,:,t) * (I - K(:,:,t) * C(:,:,t))');
 
end


%% Backward Loop: Disturbance Smoother
StT(:,T) = C(:,:,T)' * (SigmaYttm1(:,:,T) \ Ytilde(:,T)); 
XtT(:,T) = Xtt(:,T);
if nargout > 1
   SigmatT        = zeros(Nx, Nx, T);
   N              = C(:,:,T)' * (SigmaYttm1(:,:,T) \ C(:,:,T));
   SigmatT(:,:,T) = Sigmatt(:,:,T);
   % checkdiff(SigmatT(:,:,T), Sigmattm1(:,:,T) - Sigmattm1(:,:,T)  * N * Sigmattm1(:,:,T));
end

for t = T - 1 : -1 : 1
    Atilde   = A(:,:,t+1) * (I - KC(:,:,t));
    StT(:,t) = C(:,:,t)' * (SigmaYttm1(:,:,t) \ Ytilde(:,t)) + Atilde' * StT(:,t+1);
    XtT(:,t) = Xttm1(:,t) + Sigmattm1(:,:,t) * StT(:,t);
    if nargout > 1
        N               = C(:,:,t)' * (SigmaYttm1(:,:,t) \ C(:,:,t)) + Atilde' * N * Atilde;
        SigmatT(:,:,t)  = Sigmattm1(:,:,t) - Sigmattm1(:,:,t)  * N * Sigmattm1(:,:,t);
    end
end
