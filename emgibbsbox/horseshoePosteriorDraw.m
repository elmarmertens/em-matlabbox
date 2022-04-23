function [globalStatePREV, globalScalePREV, localStatePREV, localScalePREV] = ...
    horseshoePosteriorDraw(shocks2, globalStatePREV, globalScalePREV, ~, localScalePREV, dim)
% HORSESHOEPOSTERIORDRAWS ...
%
%   ...

%% VERSION INFO
% AUTHOR    : Elmar Mertens
% $DATE     : 20-Apr-2022 13:37:06 $
% $Revision : 1.00 $
% DEVELOPED : 9.12.0.1884302 (R2022a)
% FILENAME  : horseshoePosteriorDraws.m

%% prepare
if nargin < 6 || isempty(dim)
    dim = 1;
end

T       = size(shocks2, dim);
Tp1half = (T + 1) * .5;

igamrnd = @(alpha,beta) 1 ./ gamrnd(alpha, 1 ./ beta);



%% compute
this           = 1 ./ localScalePREV + .5 * shocks2 / globalStatePREV;
localStatePREV = igamrnd(1, this);
localScalePREV = igamrnd(1, 1 + 1 ./ localStatePREV);

this            = 1 ./ globalScalePREV + .5 * sum(shocks2 ./ localStatePREV, dim);
globalStatePREV = igamrnd(Tp1half, this);
globalScalePREV = igamrnd(1, 1 + 1 ./ globalStatePREV);
