function fh = newfigure(name)
% NEWFIGURE
% creates a new figure window 
% USAGE: newfigure(name) 
%        where name is an optional windowtitle 
%        (otherwise, will be named by the calling mfilename)
%
%   Coded by  Elmar Mertens, em@elmarmertens.com

this = figure;
clf reset
if nargin < 1
   m = mfilenamecaller;
   if ~isempty(m)
      set(this, 'name', m)
   end
else
   set(this, 'name', name)
end

if nargout > 0
    fh = this;
end


h = rotate3d;
set(h, 'rotatestyle', 'box');
set(this, 'defaultLegendAutoUpdate','off')

% set(gcf, 'Renderer', 'painters')