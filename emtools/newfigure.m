function fh = newfigure(name)
% NEWFIGURE
% creates a new figure window 
% USAGE: newfigure(name) 
%        where name is an optional windowtitle 
%        (otherwise, will be named by the calling mfilename)
%
% implements also the bugfix for plotting patches under unix

%   Coded by  Elmar Mertens, em@elmarmertens.com

figure
clf reset
if nargin < 1
   m = mfilenamecaller;
   if ~isempty(m)
      set(gcf, 'name', m)
   end
else
   set(gcf, 'name', name)
end

if nargout > 0
    fh = gcf;
end

% prior to rh5, we needed the following bugfix 
% ... to avoid crashes after plotting patches with large faces
% if isFRB && isunix
%    set(gcf, 'renderer', 'zbuffer')
% end

% set(gca, 'fontsize', getpref('embox', 'fontsize'))

h = rotate3d;
set(h, 'rotatestyle', 'box');
% set(gcf, 'Renderer', 'painters')