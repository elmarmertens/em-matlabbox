function [fev, fevcond] = model2fev(A,B,C,k,ndxImpulse)
% MODEL2FEV
%  computes conditional and unconditional forecast-error variances responses for state space model with matrices A,B,C
%  outputs have dimension Ny x Ny x k and Ny x Ny x k x Nw
%
% USAGE: [fev, fevcond] = model2fev(A,B,C,k)
% defaults: k = 16
%
% see also model2irf and model2acf

%   Coded by  Elmar Mertens, em@elmarmertens.com


error(abcdchk(A,B,C)) % check dimensional consistency

% nx = size(A, 1);
Nw = size(B, 2);
Ny = size(C, 1);
if nargin < 4 || isempty(k)
    k = 16;
end
if nargin < 5 || isempty(ndxImpulse)
    ndxImpulse = 1 : Nw;
end


%% unconditional FEV
fev = NaN(Ny, Ny, k);

BB  = B * B';

j = 1;
PHI = BB;
fev(:,:,j) = C * PHI * C';

for j = 2 : k
    PHI        = BB + A * PHI * A';
    fev(:,:,j) = C * PHI * C';
end

%% conditional fev
if nargout >1
    fevcond = NaN(Ny, Ny, k, length(ndxImpulse));

    for n = ndxImpulse

        BB  = B(:,n) * B(:,n)';

        j = 1;
        PHI = BB;
        fevcond(:,:,j,n) = C * PHI * C';

        for j = 2 : k
            PHI              = BB + A * PHI * A';
            fevcond(:,:,j,n) = C * PHI * C';
        end
    end
end