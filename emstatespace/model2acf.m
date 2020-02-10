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
acfy = zeros(ny, ny, lags + 1, nw);
for w = 1 : nw
    %    Js = diag(1:nw == w);
    %    acfx(:,:,1, w)   = disclyap(A, B(:,w) * B(:,w)');
    acfx(:,:,1, w)   = dlyapdoubling(A, B(:,w) * B(:,w)');
    % checkdiff(disclyap(A, B * Js * B'), dlyap(A, B * Js * B'))
    acfy(:,:,1, w)   = C * acfx(:,:,1,w) * C';
    for k = 1 : lags
        acfx(:,:,k + 1, w)   = A * acfx(:,:,k,w);
        acfy(:,:,k + 1, w)   = C * acfx(:,:,k + 1,w) * C';
    end
end

if doacc
    acfy = acf2acc(acfy);
    if nargout > 1
        acfx = acf2acc(acfx);
    end
end
        