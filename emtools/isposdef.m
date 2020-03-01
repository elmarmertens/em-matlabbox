function out = isposdef(a,semiflag)
% PURPOSE: test a matrix for Positive (semi) definiteness, using cholesky
% or -- when necessary -- using an eigenvalue decomposition
% ----------------------------------------------------------------
% USAGE: ans = isposdef(x,semiflag)
% where: x   = input matrix
%---------------------------------------------------
% RETURNS:
%        ans = 1, positive (semi) definite
%        ans = 0, not positive (semi) definite
% ----------------------------------------------------------------
% NOTES:
% if semiflag = false (default), the function tests for strict positive
% definiteness using [R,p] = chol(a); ans = (p == 0); which is fairly fast
% if semiflag = true, the function tests for positive semi-definiteness
% using an eigenvalue decomposition 
% 
% In both cases, the implementation proceeds by calling 
% [R,p] = cholcov(a,semiflag);
% ----------------------------------------------------------------
%
% See also cholcov

% Adapted from isposdef, which was written by:
% James P. LeSage, Dept of Economics
% University of Toledo
% 2801 W. Bancroft St,
% Toledo, OH 43606
% jlesage@spatial-econometrics.com

if nargin < 2 || isempty(semiflag)
    semiflag = false;
end

[R,p] = cholcov(a, semiflag); %#ok<ASGLU>
out = (p == 0);
