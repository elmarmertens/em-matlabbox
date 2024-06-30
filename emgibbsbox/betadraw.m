function b = betadraw(alpha, beta, Ndraws, rndStream)
% BETADRAW ...
%
%   ...

if nargin < 3 || isempty(Ndraws)
    Ndraws = 1;
end
if nargin < 4 || isempty(rndStream)
    rndStream = getDefaultStream;
end

%% VERSION INFO
% AUTHOR    : Elmar Mertens
% $DATE     : 12-Apr-2017 13:24:22 $
% $Revision : 1.00 $
% DEVELOPED : 9.2.0.538062 (R2017a)
% FILENAME  : betadraw.m

adraws = randn(rndStream, 2 * alpha, Ndraws);
bdraws = randn(rndStream, 2 * beta, Ndraws);

achi2 = sum(adraws.^2,1);
bchi2 = sum(bdraws.^2,1);

b     = achi2 ./ (achi2 + bchi2);

