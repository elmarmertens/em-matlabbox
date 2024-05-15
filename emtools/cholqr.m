function sqrtCC = cholqr(C, nnd)
% CHOLQR ... 
%
% sqrtCC = cholqr(C,nnd) 
% produces lower triangular sqrtCC s.t. sqrtCC * sqrtCC' = C * C'
% 
% if ndd is set to true (default) diagonals are set to non-negative values
%   ... 

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 03-May-2020 18:50:58 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 9.8.0.1359463 (R2020a) Update 1 
% FILENAME  : cholqr.m 

if nargin < 2
    nnd = true;
end


n = size(C,1);
[~, R] = qr(C');
sqrtCC = R(1:n,1:n)';

if nnd % enforce non-negative diagonals
    sqrtCC = sqrtCC * diag(sign(diag(sqrtCC)));
end
