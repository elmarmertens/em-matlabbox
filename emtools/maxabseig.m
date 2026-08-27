function rho = maxabseig(A, doDense)
% MAXABSEIG spectral radius (largest absolute eigenvalue) of a square matrix
%
% rho = maxabseig(A) returns max(abs(eig(A)))
%
% rho = maxabseig(A, doDense) with doDense = true forces the dense eig path
%
% Uses dense eig below a size threshold (~500): it is exact and, at that
% scale, faster than iterative alternatives. For larger matrices, uses
% sparse eigs('largestabs') asking for several eigenvalues, with a dense
% fallback when eigs does not converge.
%
% CAVEAT on the eigs path: the eigs convergence test certifies that the
% returned Ritz values are eigenvalues, not that they are the LARGEST ones.
% With clustered top moduli (typical for persistent VAR companion forms),
% eigs can converge to subdominant eigenvalues; a large Krylov subspace
% (SubspaceDimension) makes this much less likely but cannot rule it out.
% Pass doDense = true when correctness matters more than speed.
%
% See also eig, eigs

%   Coded by  Elmar Mertens, em@elmarmertens.com

if nargin < 2
    doDense = false;
end

n = size(A,1);

if doDense || n <= 500 % dense eig is exact and faster below the crossover
    rho = max(abs(eig(full(A))));
    return
end

try
    % the safeguard against missing the dominant eigenvalue is a LARGE
    % SUBSPACE (captures the top cluster), not a large k: demanding many
    % converged Ritz values mostly triggers non-convergence (NaN) on
    % clustered spectra; the relaxed tolerance is ample for max-root tests
    rhos = abs(eigs(sparse(A), min(5, n-1), 'largestabs', ...
        'SubspaceDimension', min(60, n), 'Tolerance', 1e-10));
    if isempty(rhos) || any(isnan(rhos))
        % unconverged eigenvalues return as NaN (with warning); any of
        % them could be the dominant one, so fall back to dense
        rho = max(abs(eig(full(A))));
    else
        rho = max(rhos);
    end
catch
    rho = max(abs(eig(full(A))));
end
