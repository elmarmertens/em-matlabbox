function [h, hanni] = plotCIbands(x, tails)
% function plotCI(x, tails)

%   Coded by  Elmar Mertens, em@elmarmertens.com

[T, N] = size(tails); %#ok<ASGLU>
if nargin < 2
	x = 1 : T;
end

x = x(:);
	
tails = sort(tails, 2); % note: this is just a crude swap of columns. it relies on tails being sortable

% if size(unique(i, 'rows'), 1) > 1
% 	error('tails sort not simply swapping columns')
% end

cla % CHECKME: really never needed?
p = plot(x, tails);
YLIM = ylim;
delete(p);

% denan
nanny = ~any(isnan(tails), 2);
tails = tails(nanny,:);
x     = x(nanny);

hold on

hanni = area(x, [tails(:,1) diff(tails, 1, 2)], min(YLIM), 'EdgeColor', 'none');

set(hanni(1), 'facecolor', ones(1,3));

switch (length(hanni) - 1) 
   case 3
      areacolors = [.8 .4 .8];
   case 7
      areacolors = [.8 .6 .4 .2 .4 .6 .8];
   case 5
      areacolors = [.8 .6 .4 .6 .8];
      % areacolors = [.75 .5 0 .5 .75];
   case 1
      areacolors = .8;
   otherwise
      error('unprepared for this number of tails ...')
end

for n = 2 : length(hanni)
   set(hanni(n), 'facecolor', repmat(areacolors(n - 1),1,3));
end

if nargout > 0
    h = hanni;
end
