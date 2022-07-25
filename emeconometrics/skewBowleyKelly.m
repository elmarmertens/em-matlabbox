function s = skewBowleyKelly(x, dim, doBowley)
% SKEWBOWLYKELLY ... 
%  
% s = skewBowleyKelly(x, dim, doBowley)
%   
% Set doBowley to true for Bowley measure (based on quartiles), 
% if otherwise uses Kelly measure (based on first and last decile and median)
%
% See also skewness

if nargin < 2 || isempty(dim)
    dim = 1;
end
if nargin < 3 || isempty(doBowley)
    doBowley = true;
end

if doBowley
    percvec = [25 50 75];
else % Kelly
    percvec = [10 50 90];
end

percentiles = prctile(x, percvec, dim);

% prepare call to subsref (to operate along arbitrary dimension dim)
s.type = '()';
s.subs = repmat({':'}, 1, ndims(x));

s.subs{dim} = 1;
Q1 = subsref(percentiles, s);
s.subs{dim} = 2;
Q2 = subsref(percentiles, s);
s.subs{dim} = 3;
Q3 = subsref(percentiles, s);

s = (Q3 + Q1 - 2 .* Q2) ./ (Q3 - Q1);

