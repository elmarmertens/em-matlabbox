function [acfy, acfx] = model2acf(A,B,C,lags,doacc)
% MODEL2ACF
% computes conditional autocovariance functions with "lags" lags for state space model with matrices A,B,C
% USAGE: function [acfy, acfx] = model2acf(A,B,C,lags)

%   Coded by  Elmar Mertens, em@elmarmertens.com

% Elmar Mertens
% www.elmarmertens.ch

if nargin < 4 || isempty(lags)
    lags = 16;
end
if nargin < 5
    doacc = false;
end

error(abcdchk(A,B,C)) % check dimensional consistency

nx = size(A, 1);
nw = size(B, 2);
ny = size(C, 1);

acfx = zeros(nx, nx, lags + 1, nw);

% lag-0 covariances per shock: joint Lyapunov solves with shared doubling
Q = zeros(nx, nx, nw);
for w = 1 : nw
    Q(:,:,w) = B(:,w) * B(:,w)';
end
acfx(:,:,1,:) = reshape(dlyapdoublingpages(A, Q), nx, nx, 1, nw);

for k = 1 : lags
    acfx(:,:,k + 1, :) = pagemtimes(A, acfx(:,:,k,:));
end
acfy = pagemtimes(pagemtimes(C, acfx), 'none', C, 'ctranspose');

if doacc
    acfy = acf2acc(acfy);
    if nargout > 1
        acfx = acf2acc(acfx);
    end
end
        