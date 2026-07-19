%% ------------------------------------------------------------------------
% function selftestBVARjeffries()
% SELFTESTBVARJEFFRIES verifies the postllf formula
%
% Test 1 (chain-rule identity, exact): for any split of the sample with
% T0 - K >= N,
%     logm(Y_{1:T}) = logm(Y_{1:T0}) + log p(Y_{T0+1:T} | Y_{1:T0})
% where the predictive is matricvariate-t:
%     p(Y1|Y0) = pi^(-N*T1/2) * GammaN((nu0+T1)/2)/GammaN(nu0/2)
%                * |V|^(-N/2) * |S0|^(nu0/2) * |S1|^(-(nu0+T1)/2)
% with nu0 = T0-K, V = I + X1*inv(X0'X0)*X1', E = Y1 - X1*PAI0,
% S1 = S0 + E'*inv(V)*E. The test computes both sides independently
% (LHS via two calls to BVARjeffries, RHS via direct linear algebra);
% each of the following bugs breaks the identity in a distinct way:
% wrong gammaln indexing, T-instead-of-dof in the pi exponent, missing
% log(2) from the IW normalization. (The GammaN constant pi^(N(N-1)/4)
% cancels in the ratio and is pinned down by Test 2 instead.)
%
% Test 2 (level check, numerical): for N = K = 1, compare postllf against
% brute-force quadrature of int int p(Y|b,s2) * s2^(-1) db ds2.

s = RandStream('mt19937ar') % , 'Seed', 42); % local stream; leaves global rng untouched

%% Test 1: chain-rule identity
T  = 60;
K  = 7;
N  = 3;
T0 = 20; % must satisfy T0 - K >= N

X = randn(s, T, K);
B = randn(s, K, N) / sqrt(K);
C = chol(0.5 * eye(N) + 0.5 * ones(N), 'lower');
Y = X * B + randn(s, T, N) * C';

[~,~,~,~, llfFull] = BVARjeffries(Y, X);
[~,~,~,~, llfSub]  = BVARjeffries(Y(1:T0,:), X(1:T0,:));

X0 = X(1:T0,:);       Y0 = Y(1:T0,:);
X1 = X(T0+1:end,:);   Y1 = Y(T0+1:end,:);
T1 = T - T0;

XX0  = X0' * X0;
PAI0 = XX0 \ (X0' * Y0);
S0   = (Y0 - X0 * PAI0)' * (Y0 - X0 * PAI0);
nu0  = T0 - K;
V    = eye(T1) + X1 * (XX0 \ X1');
E    = Y1 - X1 * PAI0;
S1   = S0 + E' * (V \ E);

logGammaRatio = sum(gammaln((nu0 + T1 + 1 - (1:N)) / 2) ...
    - gammaln((nu0 + 1 - (1:N)) / 2));
logpred = - N * T1 / 2 * log(pi) + logGammaRatio ...
    - N / 2 * logdetPD(V) + nu0 / 2 * logdetPD(S0) ...
    - (nu0 + T1) / 2 * logdetPD(S1);

err1 = llfFull - (llfSub + logpred);
fprintf('BVARjeffries self-test 1 (chain rule):  residual = %10.2e\n', err1)

%% Test 2: quadrature level check (N = K = 1)
T2 = 6;
x  = randn(s, T2, 1);
y  = x * 0.8 + 0.5 * randn(s, T2, 1);

[~,~,~,~, llf1] = BVARjeffries(y, x);

bhat = x \ y;
ssr  = sum((y - x * bhat).^2);
bse  = sqrt(ssr / (T2 - 1) / sum(x.^2));

fun  = @(b, s2) (2 * pi * s2).^(-T2/2) ...
    .* exp(- (ssr + sum(x.^2) * (b - bhat).^2) ./ (2 * s2)) ...
    .* s2.^(-1); % note: ||y-xb||^2 = ssr + x'x (b-bhat)^2

mnum = integral2(fun, bhat - 10*bse, bhat + 10*bse, ssr/500, 500*ssr, ...
    'AbsTol', 1e-14, 'RelTol', 1e-10);

err2 = llf1 - log(mnum);
fprintf('BVARjeffries self-test 2 (quadrature):  residual = %10.2e\n', err2)

%% verdict
tol1 = 1e-8;  % exact identity, machine precision expected
tol2 = 1e-3;  % numerical integration with truncated domain
if abs(err1) > tol1
    error('BVARjeffries:selftest', 'chain-rule identity violated (residual %g)', err1)
end
if abs(err2) > tol2
    error('BVARjeffries:selftest', 'quadrature level check failed (residual %g)', err2)
end
fprintf('BVARjeffries self-test: PASSED\n')


%% ------------------------------------------------------------------------
function ld = logdetPD(A)
% LOGDETPD log determinant of a symmetric positive definite matrix
ld = 2 * sum(log(diag(chol(A))));
end
