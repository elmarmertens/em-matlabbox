function [Xdraws, disturbanceDraws, X0draws, noiseDraws] = ...
    a2b2c2DisturbanceSmoothingSampler1draw(A, B, C, Ydata, X00, cholSigma00, ...
	sqrtR, rndStream)
% ABCDISTURBANCESMOOTHINGSAMPLER
% ....

% accepts only 2D inputs ABC; only sqrtR can be 3D 

%   Coded by  Elmar Mertens, em@elmarmertens.com


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


%% allocate memory
Ctilde                      = NaN(Ny,Nx,T);
[Sigmattm1, Atilde]         = deal(zeros(Nx, Nx, T));
Ztilde                      = zeros(Ny, T);
[XtT, Xttm1, Xplus]         = deal(zeros(Nx, T));

if nargout > 3 && ~isempty(sqrtR)    
    doNoiseDraws = true; 
    Ytilde       = zeros(Ny, T); 
else
    doNoiseDraws = false;
    noiseDraws   = [];
end

%% generate plus data

disturbanceplus = B * randn(rndStream, Nw, T);

if ~isempty(sqrtR) 
    Nmeasurmentnoise = size(sqrtR, 2);
    eplus            = randn(rndStream, Nmeasurmentnoise, T);
end
X0plus = X00 + cholSigma00 * randn(rndStream, Nx, 1); 

%% Forward Loop: Kalman Forecasts
[Sigma00, Sigmatt] = deal(cholSigma00 * cholSigma00');
Xtt     = zeros(Nx,1); % use zeros, since projection on difference between Y and Yplus

if isempty(sqrtR)
    noiseplus        = [];
else
    noiseplus        = zeros(Ny, T);
end

BB = B*B';

for t = 1 : T
    
    
    if t == 1
        Xplus(:,t)       = A * X0plus + disturbanceplus(:,t);
        Sigmattm1(:,:,t) = A * Sigmatt * A' + BB;
    else
        Xplus(:,t)       = A * Xplus(:,t-1) + disturbanceplus(:,t);
        Sigmattm1(:,:,t) = Atilde(:,:,t-1) * Sigmattm1(:,:,t-1) * Atilde(:,:,t-1)' + BB;
    end
    
    % priors
    Xttm1(:,t)              = A * Xtt;
    
    
    
    % observed innovation
    if isempty(sqrtR)
        
        Yplus            = C * Xplus(:,t);
        SigmaYttm1       = C * Sigmattm1(:,:,t) * C';
        
    else
        
        noiseplus(:,t)    = sqrtR(:,:,t) * eplus(:,t); %#ok<AGROW>
        Yplus             = C * Xplus(:,t) + noiseplus(:,t);
        SigmaYttm1        = C * Sigmattm1(:,:,t) * C' + sqrtR(:,:,t) * sqrtR(:,:,t)';
        
    end
    
    ytilde                  = Ydata(:,t)  - Yplus  - C * Xttm1(:,t);
    if doNoiseDraws
        Ytilde(:,t) = ytilde;
    end

    sqrtSigmaYttm1          = chol(SigmaYttm1, 'lower');
    Ztilde(:,t)             = sqrtSigmaYttm1 \ ytilde;
    Ctilde(:,:,t)           = sqrtSigmaYttm1 \ C;

    % Kalman Gain
    Ktilde                  = Sigmattm1(:,:,t) * Ctilde(:,:,t)';
    Atilde(:,:,t)           = A - A * Ktilde * Ctilde(:,:,t); % A * (I - Ktilde * Ctilde)
    
    % posteriors
    Xtt                     = Xttm1(:,t) + Ktilde * Ztilde(:,t);
   
end

%% Backward Loop: Disturbance Smoother
XtT(:,T)        = Xtt;

StT             = Ctilde(:,:,T)' * Ztilde(:,T);

if nargout > 1
    disturbancetT               = zeros(Nx, T);
    disturbancetT(:,T)          = BB * StT;
else
    disturbancetT        = [];
end

if doNoiseDraws
    noisetT            = zeros(Ny, T);
    noisetT(:,T)       = Ytilde(:,T) - C * (XtT(:,T) - Xttm1(:,T));
else
    noisetT  = [];
end


for t = (T-1) : -1 : 1
    StT         = Atilde(:,:,t)' * StT + Ctilde(:,:,t)' * Ztilde(:,t);
    XtT(:,t)    = Xttm1(:,t) + Sigmattm1(:,:,t) * StT;
    
    if ~isempty(disturbancetT)
        disturbancetT(:,t)        = BB * StT;
    end

    if doNoiseDraws
        noisetT(:,t)       = Ytilde(:,t) - C * (XtT(:,t) - Xttm1(:,t));
    end
    
end

%% sample everything together (and reorder output dimensions)
Xdraws  = Xplus + XtT;

if nargout > 1
    
    disturbanceDraws  = disturbanceplus + disturbancetT;
    
    if nargout > 2
        
        X0T      = Sigma00 * A' * StT; % note: no mean added to X0T since it is already included in X0plus
        X0draws  = X0plus + X0T;
        
        
        if doNoiseDraws % nargout > 3
            noiseDraws  = noiseplus + noisetT;
        end
    end
end
