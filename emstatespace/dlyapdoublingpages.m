function X = dlyapdoublingpages(A, Q, tolerance, maxiter)
% DLYAPDOUBLINGPAGES solves multiple discrete-time Lyapunov equations with common transition
%
% X = dlyapdoublingpages(A, Q) solves X(:,:,n) = A * X(:,:,n) * A' + Q(:,:,n)
% for every page n of the 3D array Q, via a doubling algorithm that shares
% the doubling powers of A across all pages
%
% Optional parameters: tolerance and maxiter
% The default values are tolerance = 1e-12 and maxiter = 1e4
%
% See also dlyapdoubling

%   Coded by  Elmar Mertens, em@elmarmertens.com

if nargin < 3 || isempty(tolerance)
    tolerance = 1e-12;
end
if nargin < 4 || isempty(maxiter)
    maxiter = 1e4;
end

X0   = Q;
Apow = A;

converged   = false;
iter        = 0;

while ~converged && iter < maxiter

    iter = iter + 1;

    X    = X0 + pagemtimes(pagemtimes(Apow, X0), 'none', Apow, 'ctranspose');
    Apow = Apow * Apow;     % this is the doubling step (shared across pages)

    dX          = X - X0;
    converged   = max(abs(dX(:))) < tolerance;
    X0          = X;

    if any(isnan(X0(:)))
        break
    end

end

if ~converged
    error('dlyapdoublingpages failed to converge after %d iterations.\n\t largest dX was %e\n\t tolerance %e\n\nLargest eigenvalue of transition is %e\n\n', iter, max(abs(dX(:))), tolerance, max(abs(eig(A))))
end
