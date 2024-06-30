function LLF = abcKalmanLike(A, B, C, Ydata, X00, cholSigma00, sqrtR)
% abcKalmanLike
% ....

%   Coded by  Elmar Mertens, em@elmarmertens.com

%% parse inputs
Nx                = size(A, 1);
[Ny, T]           = size(Ydata);

if nargin < 7
    sqrtR = [];
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
end

I                 = eye(Nx);


%% Forward Loop: Kalman Forecasts
Sigmatt = cholSigma00 * cholSigma00';
Xtt     = X00;
LLF     = 0;
% llf     = NaN(T,1);

for t = 1 : T

    % priors
    Sigmattm1        = A(:,:,t) * Sigmatt * A(:,:,t)' + B(:,:,t) * B(:,:,t)';
    Xttm1            = A(:,:,t) * Xtt;

    % observed innovation
    if isempty(sqrtR)
        SigmaYttm1       = C(:,:,t) * Sigmattm1 * C(:,:,t)';
    else
        SigmaYttm1        = C(:,:,t) * Sigmattm1 * C(:,:,t)' + ...
            sqrtR(:,:,t) * sqrtR(:,:,t)';
    end

    sqrtSigmaYttm1         = chol(SigmaYttm1, 'lower');
    logdetY                = 2 * sum(log(diag(sqrtSigmaYttm1)));
    Ztilde                 = sqrtSigmaYttm1 \ (Ydata(:,t)  - C(:,:,t) * Xttm1);

    LLF                    = LLF + logdetY + sum(Ztilde.^2);
    %     llf(t)                 = -.5 * (Ny * log(2 * pi) +logdetY + sum(Ztilde.^2));

    % Kalman Gain
    Ctilde                 = sqrtSigmaYttm1 \ C(:,:,t);
    Ktilde                 = Sigmattm1 * Ctilde';
    ImKC                   = I - Ktilde * Ctilde;

    %     K                      = Sigmattm1 * C(:,:,t)' / SigmaYttm1;
    %     checkdiff(K, Ktilde / sqrtSigmaYttm1);
    %     checkdiff(K * C(:,:,t), Ktilde * Ctilde);

    % posteriors
    Sigmatt                 = ImKC * Sigmattm1 * ImKC'; % Joseph form for numerical stability
    Xtt                     = Xttm1 + Ktilde * Ztilde;

end

LLF = -.5 * (Ny * log(2 * pi) * T + LLF);
% checkdiff(LLF, sum(llf));
