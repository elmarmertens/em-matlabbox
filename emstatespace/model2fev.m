function [fev, fevcond] = model2fev(A,B,C,k)
% MODEL2FEB
%  computes conditional and unconditional forecast-error variances responses for state space model with matrices A,B,C
%  outputs have dimension Ny x Ny x k and Ny x Ny x k x Nw
%
% USAGE: [fev, fevcond] = model2fev(A,B,C,k)
% defaults: k = 16
%
% see also model2irf and model2acf

%   Coded by  Elmar Mertens, em@elmarmertens.com

if nargin < 4 || isempty(k)
    k = 16;
end

error(abcdchk(A,B,C)) % check dimensional consistency

% nx = size(A, 1);
nw = size(B, 2);
ny = size(C, 1);


%% unconditional FEV
fev = NaN(ny, ny, k);

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
    fevcond = NaN(ny, ny, k, nw);

    for n = 1 : Nw
        BB  = B(:,n) * B(:,n)';

        j = 1;
        PHI = BB;
        fev(:,:,j) = C * PHI * C';

        for j = 2 : k
            PHI              = BB + A * PHI * A';
            fevcond(:,:,j,n) = C * PHI * C';
        end
    end
end