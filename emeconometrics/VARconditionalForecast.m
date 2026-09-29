function [fcstYdraws, zdraws, nushocks, fcstYhat, info] = VARconditionalForecast(PAI, Xjumpoff, sqrtInnov, Ycond, Ndraws, rndStream, opts)
% VARCONDITIONALFORECAST Draw VAR forecasts conditional on observable paths.
%
%   [fcstYdraws, zdraws, nushocks, fcstYhat, info] = VARconditionalForecast( ...
%       PAI, Xjumpoff, sqrtInnov, Ycond, Ndraws, rndStream, opts)
%
%   This function implements hard conditional forecasting with all future
%   shocks available to satisfy the conditions. It does not identify or
%   restrict a subset of structural shocks.
%
%   Model notation follows the CCMM VAR samplers:
%
%       X(t+h) = fcstA * X(t+h-1) + fcstB * sqrtInnov(:,:,h) * z(t+h)
%       Y(t+h) = X(t+h, ndfcstY)
%
%   where PAI is K x N, Xjumpoff is K x 1 ordered as
%   [constant; Y(t); Y(t-1); ...; Y(t-p+1)], and z(t+h) ~ N(0,I).
%
%   Inputs
%   ------
%   PAI        : K x N VAR coefficient draw, ordered as in the CCMM code.
%   Xjumpoff   : K x 1 forecast jump-off state.
%   sqrtInnov  : N x N or N x N x H matrix mapping standard shocks into
%                reduced-form innovations. For the constant-VCV model use
%                cholSigma. For SV models pass invA_ * diag(futureSV(:,h)).
%   Ycond      : N x H matrix of conditioning values. Use NaN for entries
%                that are not conditioned upon.
%   Ndraws     : number of conditional forecast draws.
%   rndStream  : optional RandStream.
%   opts       : optional struct with fields:
%                .pinvTolerance
%                .constraintTolerance
%                .strict
%
%   Outputs
%   -------
%   fcstYdraws : N x H x Ndraws conditional forecast draws.
%   zdraws     : N x H x Ndraws standardized future shocks.
%   nushocks   : N x H x Ndraws reduced-form future innovations.
%   fcstYhat   : N x H conditional mean path.
%   info       : diagnostics and conditioning details.

if nargin < 5 || isempty(Ndraws)
    Ndraws = 1;
end
if nargin < 6 || isempty(rndStream)
    if exist('getDefaultStream', 'file') == 2
        rndStream = getDefaultStream;
    else
        rndStream = RandStream.getGlobalStream;
    end
end
if nargin < 7 || isempty(opts)
    opts = struct();
end

if ~isfield(opts, 'pinvTolerance') || isempty(opts.pinvTolerance)
    opts.pinvTolerance = 1e-10;
end
if ~isfield(opts, 'constraintTolerance') || isempty(opts.constraintTolerance)
    opts.constraintTolerance = 1e-8;
end
if ~isfield(opts, 'strict') || isempty(opts.strict)
    opts.strict = true;
end

[K, N] = size(PAI);
if size(Xjumpoff, 1) ~= K || size(Xjumpoff, 2) ~= 1
    error('Xjumpoff must be K x 1, with K = size(PAI,1).');
end
if size(Ycond, 1) ~= N
    error('Ycond must be N x H, matching the CCMM forecast-array orientation.');
end

H = size(Ycond, 2);
p = (K - 1) / N;
if abs(p - round(p)) > eps
    error('The number of rows in PAI must equal 1 + N * p.');
end
p = round(p);

sqrtInnov = localExpandSqrtInnov(sqrtInnov, N, H);

[fcstA, fcstB, ndfcstY] = localCompanionForm(PAI, N, p);
[muY, impactY] = localForecastMoments(fcstA, fcstB, ndfcstY, Xjumpoff, sqrtInnov, N, H);

muStack = muY(:);
yCondStack = Ycond(:);
conditionMask = isfinite(yCondStack);
numConstraints = sum(conditionMask);

NH = N * H;
if numConstraints == 0
    zStackDraws = randn(rndStream, NH, Ndraws);
    zStackMean = zeros(NH, 1);
    nullBasis = eye(NH);
    residualMean = [];
else
    D = impactY(conditionMask, :);
    r = yCondStack(conditionMask) - muStack(conditionMask);

    DD = D * D';
    DD = (DD + DD') * 0.5;
    zStackMean = D' * (pinv(DD, opts.pinvTolerance) * r);
    residualMean = D * zStackMean - r;

    if max(abs(residualMean)) > opts.constraintTolerance
        msg = sprintf('Conditional path is not exactly feasible; max residual is %.3g.', max(abs(residualMean)));
        if opts.strict
            error('%s', msg);
        else
            warning('%s', msg);
        end
    end

    nullBasis = null(D);
    zStackDraws = repmat(zStackMean, 1, Ndraws);
    if ~isempty(nullBasis)
        zStackDraws = zStackDraws + nullBasis * randn(rndStream, size(nullBasis, 2), Ndraws);
    end
end

fcstYdraws = NaN(N, H, Ndraws);
zdraws = NaN(N, H, Ndraws);
nushocks = NaN(N, H, Ndraws);

YstackDraws = repmat(muStack, 1, Ndraws) + impactY * zStackDraws;
for ndraw = 1:Ndraws
    fcstYdraws(:, :, ndraw) = reshape(YstackDraws(:, ndraw), N, H);
    zdraws(:, :, ndraw) = reshape(zStackDraws(:, ndraw), N, H);
    for h = 1:H
        nushocks(:, h, ndraw) = sqrtInnov(:, :, h) * zdraws(:, h, ndraw);
    end
end

fcstYhat = reshape(muStack + impactY * zStackMean, N, H);

if numConstraints == 0
    constraintResidualByDraw = [];
else
    constraintResidualByDraw = YstackDraws(conditionMask, :) - yCondStack(conditionMask);
end

info = struct();
info.forecastType = 'conditional_observables_all_shocks';
info.N = N;
info.H = H;
info.p = p;
info.numConstraints = numConstraints;
info.conditionMask = reshape(conditionMask, N, H);
info.Ycond = Ycond;
info.fcstYuncondMean = muY;
info.fcstYcondMean = fcstYhat;
info.zMean = reshape(zStackMean, N, H);
info.nullRank = size(nullBasis, 2);
info.constraintResidualMean = residualMean;
info.constraintResidualByDraw = constraintResidualByDraw;
info.maxAbsConstraintResidual = localMaxAbs(constraintResidualByDraw);
info.maxAbsZByDraw = max(reshape(abs(zdraws), [], Ndraws), [], 1);
info.maxAbsZ = max(info.maxAbsZByDraw);
info.maxAbsNuByDraw = max(reshape(abs(nushocks), [], Ndraws), [], 1);
info.maxAbsNu = max(info.maxAbsNuByDraw);

end

function sqrtInnov = localExpandSqrtInnov(sqrtInnov, N, H)

if isequal(size(sqrtInnov), [N, N])
    sqrtInnov = repmat(sqrtInnov, 1, 1, H);
elseif ndims(sqrtInnov) == 3 && size(sqrtInnov, 1) == N && size(sqrtInnov, 2) == N && size(sqrtInnov, 3) == H
    % already in horizon-specific form
else
    error('sqrtInnov must be N x N or N x N x H.');
end

end

function [fcstA, fcstB, ndfcstY] = localCompanionForm(PAI, N, p)

K = size(PAI, 1);
fcstA = zeros(K, K);
fcstA(1, 1) = 1;

ndfcstY = 1 + (1:N);
fcstA(ndfcstY, :) = PAI';

if p > 1
    fcstA(1 + N + 1:end, 2:end) = [eye(N * (p - 1)), zeros(N * (p - 1), N)];
end

fcstB = zeros(K, N);
fcstB(ndfcstY, :) = eye(N);

end

function [muY, impactY] = localForecastMoments(fcstA, fcstB, ndfcstY, Xjumpoff, sqrtInnov, N, H)

K = size(fcstA, 1);
selectY = zeros(N, K);
selectY(:, ndfcstY) = eye(N);

muY = NaN(N, H);
impactY = zeros(N * H, N * H);

X = Xjumpoff;
for h = 1:H
    X = fcstA * X;
    muY(:, h) = X(ndfcstY);

    rowNdx = (h - 1) * N + (1:N);
    for shockH = 1:h
        colNdx = (shockH - 1) * N + (1:N);
        impactY(rowNdx, colNdx) = selectY * (fcstA ^ (h - shockH)) * fcstB * sqrtInnov(:, :, shockH);
    end
end

end

function value = localMaxAbs(x)

if isempty(x)
    value = [];
else
    value = max(abs(x(:)));
end

end
