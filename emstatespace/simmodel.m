function [Y, X] = simmodel(A, B, C, w, X0)
% SIMMODEL
% simulates ABC state space model
%
% USAGE: [Y, X] = simmodel(A, B, C, w, X0) 
% OR:    [Y, X] = simmodel(A, B, C, T, X0) 
% where A, B, C are state space matrices, T is length of sample (or: w is sequence of shocks) 
%       and X0 specifies initial conditions

error(abcdchk(A,B,C)) % check dimensional consistency

nx = size(A, 1);
nw = size(B, 2);
ny = size(C, 1); %#ok<NASGU>

if isscalar(w)
   T = w;
   w = randn(T, nw);
else
   [T, check] = size(w);
   if check ~= nw
      error('dimension of w not conistent with number of shocks')
   end
end
   
if nargin < 5
   X0 = zeros(1, nx);
else
   X0 = X0(:)'; % make sure it is a row vector, w/o dimension check
end

%% fast version
X = ltitr(A, B, [w; zeros(1, nw)], X0); % augmenting the shocks with a row of zeros is a peculiarity of ltitr, the additional row is omitted in the computation
X = X(2:end, :);
Y = X * C';
 
% %% slow and pedestrian version
% x = zeros(T, nx);
% y = zeros(T, ny);
% 
% x(1,:) = X0 * A' + w(1,:) * B';
% y(1,:) = x(1,:) * C';
% 
% for t = 2 : T
%    x(t,:) = x(t-1,:) * A' + w(t, :) * B'; % notice the transpose!
%    y(t,:) = x(t,:) * C';   
% end
% 
% 
% checkdiff(X, x);
% checkdiff(Y, y);
