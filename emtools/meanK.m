function mk = meanK(x, k, dim)
% trailing moverage average 
% USAGE: function mk = meanK(x, k, dim)
% wraps call to matlab's movmean
%

%   Coded by  Elmar Mertens, em@elmarmertens.com

narginchk(2,3)
if nargin < 3
   dim = 1;
end

mk0 = movmean(x, [k 0], dim, 'omitnan', 'endpoints', 'discard');

% pad with NaN

% prepare call to subsref (to operate along arbitrary dimension dim)
s.type = '()';
s.subs = repmat({':'}, 1, ndims(x));
s.subs{dim} = k+1 : size(x,dim);

mk = NaN(size(x));
mk = subsasgn(mk,s,mk0);
