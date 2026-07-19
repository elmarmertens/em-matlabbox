function iwishchols = iwishcholdrawbartlett(cholSigmaT, dof, Ndraw, rndStream)
% IWISHCHOLDRAWBARTLETT produces random number from inverse wishart distribution via Bartlett decomposition
%
% USAGE: iwishchols = iwishcholdrawbartlett(cholSigmaT, dof, Ndraw, rndStream)
%
% return left square root factor of IWISH (though not lower triangular)
%
% Distributionally identical to iwishcholdraw, but constructs the Cholesky
% factor of the underlying Wishart(I,dof) draw directly via the Bartlett
% decomposition: W = L * L' with L lower triangular,
%   L(i,i) = sqrt(chi2(dof-i+1)) and L(i,j) ~ N(0,1) for i > j.
% This needs only N*(N+1)/2 random numbers per draw and avoids forming the
% N x dof normals, the Gram product z*z' and its Cholesky factorization.
%
% Note: consumes a different random-number sequence than iwishcholdraw, so
% seeded results differ (draw by draw) between the two functions.
%
% Notice, for a proper draw (with first moment defined), need dof > N + 1.
%
% See also iwishcholdraw, iwishdraw

%   Coded by  Elmar Mertens, em@elmarmertens.com

%% parse inputs
if nargin < 2
    dof = [];
end
if nargin < 3
    Ndraw = 1;
end
if nargin < 4
    rndStream = getDefaultStream;
end

%% check dof
N = size(cholSigmaT, 1);
if isempty(dof)
    dof = N;
end
if dof < N
    error('em:iwishcholdrawbartlett', 'need dof >= N (dof=%d, N=%d)', dof, N)
end

%% draw random numbers
zlower   = randn(rndStream, N, N, Ndraw);                             % only strict lower triangle used
halfdofs = (dof - (0 : N-1)') / 2;
% randg does not accept a RandStream argument; swap rndStream in as global stream
prevStream = RandStream.setGlobalStream(rndStream);
chi2diag   = 2 * randg(repmat(halfdofs, 1, Ndraw));                   % chi2(dof-i+1) = 2 * gamma((dof-i+1)/2, 2)
RandStream.setGlobalStream(prevStream);

%% loop over Ndraw and construct inverse-wisharts
iwishchols = NaN(N,N,Ndraw);
for n = 1 : Ndraw
    L                 = tril(zlower(:,:,n), -1);
    L(1:N+1:end)      = sqrt(chi2diag(:,n));
    iwishchols(:,:,n) = cholSigmaT / L'; % L' is the upper choleski factor of a Wishart(I,dof) draw
end
