function sqrtCC = cholqr(C)
% CHOLQR ... 
%
% sqrtCC = cholqr(C) 
% produces lower triangular sqrtCC s.t. sqrtCC * sqrtCC' = C * C'
%  
%   ... 

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 03-May-2020 18:50:58 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 9.8.0.1359463 (R2020a) Update 1 
% FILENAME  : cholqr.m 



n = size(C,1);
[~, R] = qr(C');
sqrtCC = R(1:n,1:n)';
% enforce non-negative diagonals
sqrtCC = sqrtCC * diag(sign(diag(sqrtCC)));
