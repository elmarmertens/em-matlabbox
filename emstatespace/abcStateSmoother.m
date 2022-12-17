function [XtT, SigmatT, Xtt, Sigmatt] = abcStateSmoother(A, B, C, Ydata, X00, Sigma00, R)
% ABCSTATESMOOTHER produces smoothed estimates of state vector
% uses Carter-Kohn "state" smoothing (as in Hamilton), which is computationally ineffcient
%
% USAGE [XtT, SigmatT] = abcStateSmoother(A, B, C, Ydata, X00, Sigma00, R)
%  
% See also abcDisturbanceSmoother

%   Coded by  Elmar Mertens, em@elmarmertens.com

%% VERSION INFO
% AUTHOR    : Elmar Mertens
% FILENAME  : abcStateSmoother.m

%% parse inputs
if nargin < 7
    R = [];
end

%% init Variables and allocate memory
Nx                = size(A, 1);
[~, T]           = size(Ydata);
% Nw                = size(B,2);

Ydata             = permute(Ydata, [1 3 2]);

if ismatrix(A)
    A = repmat(A, [1 1 T]);
end
if ismatrix(B)
    B = repmat(B, [1 1 T]);
end
if ismatrix(C)
    C = repmat(C, [1 1 T]);
end
if ~isempty(R) && ismatrix(R)
    R = repmat(R, [1 1 T]);
end

Ix          = eye(Nx);

[Sigmatt, Sigmattm1]    = deal(zeros(Nx, Nx, T));
Xtt                     = zeros(Nx, T);

%% Forward Loop: Kalman Forecasts
for t = 1 : T
    
    % setup priors
    if t == 1
        Sigmattm1(:,:,t)    = A(:,:,t) *  Sigma00 * A(:,:,t)' + B(:,:,t) * B(:,:,t)';
        Xttm1               = A(:,:,t) *  X00;
    else
        Sigmattm1(:,:,t)    = A(:,:,t) *  Sigmatt(:,:,t-1) * A(:,:,t)' + B(:,:,t) * B(:,:,t)';
        Xttm1               = A(:,:,t) *  Xtt(:,t-1);
    end
    
    % observed innovation
    Ytilde                  = Ydata(:,:,t) - C(:,:,t) * Xttm1;
    
    
    % Kalman Gain
    SigmaYttm1              = C(:,:,t) * Sigmattm1(:,:,t) * C(:,:,t)';
    if ~isempty(R)        
        SigmaYttm1   = SigmaYttm1 + R(:,:,t);
    end
    K                      = (Sigmattm1(:,:,t) * C(:,:,t)') / SigmaYttm1;
    
    % posteriors
    Sigmatt(:,:,t)         = (Ix - K * C(:,:,t)) * Sigmattm1(:,:,t);
    Xtt(:,t)               = Xttm1 + K * Ytilde;
    
end


%% Kalman Backcast and State Draws
XtT       = zeros(Nx, T);
XtT(:,T)  = Xtt(:,T);

if nargout > 1
    SigmatT        = zeros(Nx, Nx, T);
    SigmatT(:,:,T) = Sigmatt(:,:,T);
else
    SigmatT = [];
end

for t = T-1 : -1 : 1
 
   Sigmatp1t      = A(:,:,t+1) * Sigmatt(:,:,t) * A(:,:,t+1)' + B(:,:,t+1) * B(:,:,t+1)';
   J              = Sigmatt(:,:,t) * A(:,:,t+1)' / Sigmatp1t;
   XtT(:,t)       = Xtt(:,t) + J * (XtT(:,t+1) - A(:,:,t+1) * Xtt(:,t));

   if ~isempty(SigmatT)
       SigmatT(:,:,t) = Sigmatt(:,:,t) + J * (SigmatT(:,:,t+1) - Sigmattm1(:,:,t+1)) * J';
   end
   
end
