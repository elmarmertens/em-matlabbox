function [Sigma, K, P, SigmaStar] = abckalman(A,B,C,useDare)
% ABCKALMAN solves a Kalman Filter equation by using matlab's dare (or VFI)
% [Sigma, K, P, SigmaStar] = abckalman(A,B,C)
% 
% where Sigma is innovation-VCV, K is gain, P is projection matrix 
% and SigmaStar is residual-VCV
% 
% See also abckalmani, abckalmanvfi

%   Coded by  Elmar Mertens, em@elmarmertens.com

%% VERSION INFO
% AUTHOR    : Elmar Mertens
% $DATE     : 04-May-2010 11:45:27 $
% $Revision : 1.00 $
% DEVELOPED : 7.8.0.347 (R2009a)
% FILENAME  : abckalman.m

if nargin < 4
    useDare = false;
end
if useDare
    try
        Sigma = dare(A', C', B*B', zeros(size(C, 1)));
    catch dareME
        warning('DARE failed, using riccati via VFI: %s', dareME.identifier)
        Sigma = riccati(A', C', B*B', zeros(size(C, 1)));
    end
else
    Sigma = riccati(A', C', B*B', zeros(size(C, 1)));
end
varz = C * Sigma * C';
covxz = Sigma * C';
if rcond(varz) > 1e-10
   K = covxz / varz;
else
   K = covxz * pinv(varz);
end

if nargout > 2
   P = K * C;
   if nargout > 3
      SigmaStar = (eye(size(A)) - P) * Sigma; % note: equal to Sigma * (eye(size(A)) - P)'
   end
end
