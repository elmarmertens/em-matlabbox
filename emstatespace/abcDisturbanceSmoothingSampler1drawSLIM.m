function [Xdraws, X0draws] = ...
    abcDisturbanceSmoothingSampler1drawSLIM(A, B, C, Ydata, X00, cholSigma00, ...
    sqrtR, rndStream)
% ABCDISTURBANCESMOOTHINGSAMPLER
% ....

%   Coded by  Elmar Mertens, em@elmarmertens.com

%% VERSION INFO
% AUTHOR    : Elmar Mertens
% $DATE     : 08-Aug-2009 17:58:16 $
% $Revision : 1.00 $
% DEVELOPED : 7.7.0.471 (R2008b)
% FILENAME  : abcDisturbanceSmoothingSampler.m


%% parse inputs
Nx                = size(A, 1);
[Ny, T]           = size(Ydata);
Nw                = size(B,2);

if nargin < 7
    sqrtR = [];
end

if nargin < 8
    rndStream = getDefaultStream;
end

%% init Variables and allocate memory
if ~isempty(sqrtR) && ismatrix(sqrtR)
    sqrtR = repmat(sqrtR, [1 1 T]);
end


if ismatrix(A)
    A = repmat(A, [1 1 T]);
end
if ismatrix(B)
    B = repmat(B, [1 1 T]);
end
if ismatrix(C)
    C = repmat(C, [1 1 T]);
    % warning('em:msg', 'C should be three dimensional when there are missing data')
end

I                 = eye(Nx);

%% allocate memory
Ctilde                      = NaN(Ny,Nx,T);
[Sigmattm1, ImKC]           = deal(zeros(Nx, Nx, T));
Ztilde                      = zeros(Ny, T);
[XtT, Xttm1, Xplus]         = deal(zeros(Nx, T));


%% generate plus data

wplus   = randn(rndStream, Nw, T);
if ~isempty(sqrtR)
    Nmeasurmentnoise = size(sqrtR, 2);
    eplus            = randn(rndStream, Nmeasurmentnoise, T);
end
X0plus = X00 + cholSigma00 * randn(rndStream, Nx, 1);

%% Forward Loop: Kalman Forecasts
[Sigma00, Sigmatt] = deal(cholSigma00 * cholSigma00');
Xtt     = zeros(Nx,1); % use zeros, since projection on difference between Y and Yplus


for t = 1 : T

    if t == 1
        Xplus(:,t) = A(:,:,t) * X0plus + B(:,:,t) * wplus(:,t);
    else
        Xplus(:,t) = A(:,:,t) * Xplus(:,t-1) + B(:,:,t) * wplus(:,t);
    end

    % priors
    Sigmattm1(:,:,t)        = A(:,:,t) * Sigmatt * A(:,:,t)' + B(:,:,t) * B(:,:,t)';
    Xttm1(:,t)              = A(:,:,t) * Xtt;



    % observed innovation
    if isempty(sqrtR)

        Yplus            = C(:,:,t) * Xplus(:,t);
        SigmaYttm1       = C(:,:,t) * Sigmattm1(:,:,t) * C(:,:,t)';

    else

        Yplus             = C(:,:,t) * Xplus(:,t) + sqrtR(:,:,t) * eplus(:,t);
        SigmaYttm1        = C(:,:,t) * Sigmattm1(:,:,t) * C(:,:,t)' + ...
            sqrtR(:,:,t) * sqrtR(:,:,t)';

    end
    sqrtSigmaYttm1          = chol(SigmaYttm1, 'lower');

    Ztilde(:,t)             = sqrtSigmaYttm1 \ (Ydata(:,t)  - Yplus  - C(:,:,t) * Xttm1(:,t));
    Ctilde(:,:,t)           = sqrtSigmaYttm1 \ C(:,:,t);

    % Kalman Gain
    Ktilde                  = Sigmattm1(:,:,t) * Ctilde(:,:,t)';
    ImKC(:,:,t)             = I - Ktilde * Ctilde(:,:,t);

    % posteriors
    Sigmatt                 = ImKC(:,:,t) * Sigmattm1(:,:,t);

    Xtt                     = Xttm1(:,t) + Ktilde * Ztilde(:,t);

end

%% Backward Loop: Disturbance Smoother
XtT(:,T)        = Xtt;

StT             = Ctilde(:,:,T)' * Ztilde(:,T);


for t = (T-1) : -1 : 1
    Atilde      = A(:,:,t+1) * ImKC(:,:,t);
    StT         = Atilde' * StT + Ctilde(:,:,t)' * Ztilde(:,t);
    XtT(:,t)    = Xttm1(:,t) + Sigmattm1(:,:,t) * StT;
end

%% sample everything together (and reorder output dimensions)
Xdraws  = Xplus + XtT;


if nargout > 1
    X0T      = Sigma00 * A(:,:,1)' * StT; % note: no mean added to X0T since it is already included in X0plus
    X0draws  = X0plus + X0T;
end
