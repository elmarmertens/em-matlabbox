function iwishchols = iwishcholdraw(cholSigmaT, dof, Ndraw, rndStream)
% IWISHCHOLDRAW produces random number from inverse wishart distribution
%
% USAGE: ...
%
% return left square root factor of IWISH (though not lower triangular)
%
% Notice, for a proper draw (with first moment defined) , need dof > N + 1.
%   ...

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

%% draw random normal numbers
z       = randn(rndStream, N, dof, Ndraw);
% todo: this could be improved w/bartlett decomposition ...

%% loop over Ndraw and construct inverse-wisharts
iwishchols = NaN(N,N,Ndraw);
for n = 1 : Ndraw
    sqrtZZ              = chol(z(:,:,n) * z(:,:,n)'); % note: absensce of transpose, using right choleski
    %     [~, sqrtZZ]         = qr(z(:,:,n)', 0); % Choleski is much quicker than QR
    iwishchols(:,:,n)   = cholSigmaT / sqrtZZ; % note: just a square root factor, not a lower triangular
    %     checkdiff(iwishchols(:,:,n) * iwishchols(:,:,n)', cholSigmaT / (z(:,:,n) * z(:,:,n)') * cholSigmaT');
end

