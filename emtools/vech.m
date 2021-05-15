function [v, ndx] = vech(m, ndx)
% VECH converts symmetric matrix into vector of uniqe elements (lower triangular)
% USAGE: [v, ndx] = vech(m, ndx)
% ndx is optional parameter of logical indices into lower diagonal elements of m
% 
% See also ivech, nvech

%   Coded by  Elmar Mertens, em@elmarmertens.com

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 28-Mar-2009 14:29:12 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 7.7.0.471 (R2008b) 
% FILENAME  : vech.m 

if nargin < 2
   ndx = logical(tril(ones(size(m))));
end
v   = m(ndx);
