function [XtT, SigmatT, xshocktT, llf, Sigmattm1] = a2b2c2noiseDisturbanceSmoother(A, B, C, Ydata, X00, Sigma00, R)
% ABCDISTURBANCESMOOTHER computes smoothed state estimates for ABC system 
% uses disturbance smoothing from Koopman & Co
%
% [XtT, SigmatT, xshocktT, llf, Sigmattm1] = abcDisturbanceSmoother(A, B, C, Ydata, X00, Sigma00, R)
% where R is optional VCV matrix of noise in observer system (default is zero)
%
% See also abcStateSmoother

%   Coded by  Elmar Mertens, em@elmarmertens.com


%% init Variables and allocate memory
Nx                = size(B,1);
% Nw                = size(B,2); 
[Ny, T, Ndraws]   = size(Ydata);


if nargin < 7
   R = [];
elseif ismatrix(R)
   R = repmat(R, [1 1 T]);
end

Ydata             = permute(Ydata, [1 3 2]);

Ix                = eye(Nx);
Iy                = eye(Ny);

[Sigmattm1, KC]             = deal(zeros(Nx, Nx, T));
invSigmaYttm1               = zeros(Ny, Ny, T);
Ytilde                      = zeros(Ny, Ndraws, T);
[XtT, Xttm1]                = deal(zeros(Nx, Ndraws, T));

BB = B * B';

if nargout > 3
   llf                      = zeros(T, Ndraws);
else 
   llf = [];
end

%% Forward Loop: Kalman Forecasts
for t = 1 : T
    
    
    % setup priors
    if t == 1
        Sigmattm1(:,:,1)    = A *  Sigma00 * A' + BB;
        Xttm1(:,:,1)        = A * repmat(X00, 1, Ndraws);
    else
        Sigmattm1(:,:,t)    = A *  Sigmatt * A' + BB;
        Xttm1(:,:,t)        = A * Xtt;
    end
    
    % observed innovation
    Ytilde(:,:,t)          = Ydata(:,:,t) - C * Xttm1(:,:,t); % note: no need for bsxfun, since Ndraw is defined by dimensions of Y
    
    % Kalman Gain
    if isempty(R)
        SigmaYttm1         = C * Sigmattm1(:,:,t) * C';
    else
        SigmaYttm1         = C * Sigmattm1(:,:,t) * C' + R(:,:,t);
    end
    
    if det(SigmaYttm1) == 0  % || rcond(SigmaYttm1) < 1e-8 % det is slightly faster ...
       invSigmaYttm1(:,:,t)    = pinv(SigmaYttm1);
    else
       invSigmaYttm1(:,:,t)    = Iy / SigmaYttm1;
    end
    % notice: any(all(C == 0, 2)); would also catch missing values,
    % but performed slowest during profiling


    K                       = (Sigmattm1(:,:,t) * C') * invSigmaYttm1(:,:,t);
    KC(:,:,t)               = K * C;
    % posteriors
    Sigmatt                 = (Ix - KC(:,:,t)) * Sigmattm1(:,:,t);
    
    %     checkdiff(Sigmatt, ...
    %         (Ix - KC(:,:,t)) * Sigmattm1(:,:,t) * (Ix - KC(:,:,t))' + K * R(:,:,t) * K', ...
    %         1e-8);
    
    Xtt                     = Xttm1(:,:,t) + K * Ytilde(:,:,t);
 
    % likelihood
    if ~isempty(llf)
       for n = 1 : Ndraws
          llf(t,n)            = -.5 * (Ny * log(2 * pi) + log(det(SigmaYttm1)) + Ytilde(:,n,t)' * invSigmaYttm1(:,:,t) * Ytilde(:,n,t));
       end
    end
    
end


%% Backward Loop: Disturbance Smoother
XtT(:,:,T)  = Xtt;
StT         = C' * (invSigmaYttm1(:,:,T) * Ytilde(:,:,T)); 

% SigmatT = NaN;
if nargout > 1
    SigmatT         = zeros(Nx, Nx, T);
    SigmatT(:,:,T)  = Sigmatt;
    
    N               = C' * invSigmaYttm1(:,:,T) * C;
    checkdiff(SigmatT(:,:,T), Sigmattm1(:,:,T) - Sigmattm1(:,:,T)  * N * Sigmattm1(:,:,T));
    
else
    SigmatT = [];
end


if nargout > 2
    xshocktT         = NaN(Nx, Ndraws, T);
    xshocktT(:,:,T)  = BB' * StT;
else
    xshocktT = [];
end

for t = T-1 : -1 : 1
    Atilde      = A * (Ix - KC(:,:,t));
    StT         = C' * (invSigmaYttm1(:,:,t) * Ytilde(:,:,t)) + ...
        Atilde' * StT;
    XtT(:,:,t)  = Xttm1(:,:,t) + Sigmattm1(:,:,t) * StT;
    
    if ~isempty(SigmatT)
         N               = C' * invSigmaYttm1(:,:,t) * C + Atilde' * N * Atilde;
         SigmatT(:,:,t)  = Sigmattm1(:,:,t) - Sigmattm1(:,:,t)  * N * Sigmattm1(:,:,t);
    end
     
    if ~isempty(xshocktT)
        xshocktT(:,:,t) = BB' * StT;
    end
end

%% reshape output
XtT = permute(XtT, [1 3 2]);
if ~isempty(xshocktT)
    xshocktT = permute(xshocktT, [1 3 2]);
end
