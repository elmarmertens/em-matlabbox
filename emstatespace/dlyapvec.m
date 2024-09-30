function X = dlyapvec(A,B,Q)
% DLYAPDOUBLING solves the discrete-time Lyapunov function with a doubling algorithm
% this code mimicks Matlab's own dlyap *for pure Lyapunov equations*
% 
% X = dlyapvec(A,Q) solves X = A * X * A' + Q
%
% X = dlyapvec(A,B,Q) solves the Sylvester equation X = A * X * B' + Q 
%
% Optional parameters: tolerance and maxiter
% which can be specified via dlyap(A,Q,tolerance) and dlyap(A,Q,tol,maxiter) 
% (or dlyap(A,B,Q,tolerance) and dlyap(A,B,Q,tolerance,maxiter)). 
% The default values are tolerance = 1e-12 and maxiter = 1e4
%
% See also dlyapdoubling

%   Coded by  Elmar Mertens, em@elmarmertens.com

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 24-Nov-2010 11:58:24 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 7.7.0.471 (R2008b) 
% FILENAME  : dlyap.m 

if nargin == 2 || isempty(Q)
    Q = B;
    B = A;
end


n           = size(Q, 1);
In2         = eye(n * n);

vecX        = (In2 - kron(B, A)) \ Q(:);
X           = reshape(vecX, n, n);

