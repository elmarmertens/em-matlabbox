function [bdraw, resid] = bayesRegressionSlopesGibbsDrawInvPrior(y, X, h, iV0b0, iV0, Ndraws, rndStream)
% BAYESREGRESSIONSLOPESGIBBSDRAW performs Gibbs steps for linear regression model 
%   this function draws slopes conditional on precision (inverse residual variance)
%   Usage: [bdraw, resid] =  bayesRegressionSlopesGibbsDrawInvPrior(y, X, h, iV0b0, iV0, Ndraws, rndStream)
%
% See also bayesRegressionGibbsDraw 

%   Coded by  Elmar Mertens, em@elmarmertens.com

%% VERSION INFO
% AUTHOR    : Elmar Mertens
% $DATE     : 17-Jan-2009 16:30:05 $
% $Revision : 1.00 $
% DEVELOPED : 7.7.0.471 (R2008b)
% FILENAME  : bayesRegressionGibbsDraw.m

if nargin < 6 || isempty(Ndraws)
   Ndraws = 1;
end
if nargin < 7 || isempty(rndStream)
   rndStream = getDefaultStream;
end

%% Some prelim preparations
[T, K]      = size(X); %#ok<ASGLU>

if nargin < 4 || isempty(iV0b0)
    iV0b0 = zeros(K, 1);
end
if nargin < 5 || isempty(iV0) % diffuse prior
    iV0 = zeros(K);
elseif isvector(iV0)
    iV0 = diag(iV0);
end
if nargin < 6 || isempty(Ndraws)
   Ndraws = 1;
end

XX    = X' * X;
Xy    = X' * y;

% posterior precision
iV     = iV0 + XX * h; 
sqrtiV = chol(iV); % upper triangular

%% draw regression slope
zdraw  = randn(rndStream, K, Ndraws);

bdraw     = sqrtiV' \ (iV0b0 + Xy * h); % preparatory step
bdraw     = sqrtiV \ (bdraw + zdraw); 

if nargout > 1
    resid       = y - X * bdraw;
end


   
