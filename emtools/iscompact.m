function flag = iscompact(nanny)
% function flag = iscompact(nanny);
% returns true if the *logical* sample index nanny is compact 
% (i.e. prunes only observations at beginning and/or end)
% flag  = ~all(diff(find(~nanny(:))) == 1)

%   Coded by  Elmar Mertens, em@elmarmertens.com

if ~islogical(nanny)
   nanny = isnan(nanny);
end
if ~isvector(nanny)
   error('nanny should be vector')
end

flag  = all(diff(find(~nanny(:))) == 1);
