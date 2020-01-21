function h = plotdraws(x, label, muflag, farbe, newfigflag)
% HISTDRAWS ... 
%  
%       h = histdraws(x, label)
%   ... 

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 12-Mar-2012 14:55:27 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 7.11.0.584 (R2010b) 
% FILENAME  : histdraws.m 

if nargin < 2 || isempty(label)
    label = '';
end

if nargin < 3 || isempty(muflag)
    muflag = true;
end

if nargin < 4 || isempty(farbe)
    farbe = [0 0 1]; % blue
end

if nargin < 5 || isempty(newfigflag)
    newfigflag = true;
end

x      = x(:);
% Ndraws = length(x);

[pdf, xtick] = ksdensity(x);

if newfigflag
    newfigure(label)
end
hanni = NaN(4,1);
hold on
plot(xtick, pdf, '-', 'linewidth', 2, 'color', farbe)
if muflag
    mu  = mean(x);
    %      mid = median(x);
    % vol = std(x);
    hanni(1) = plotvertline(mu, [], '-.', 'linewidth', 2, 'color', farbe);
    %     hanni(2) = plotvertline(mid, [], '--', 'linewidth', 2, 'color', farbe);
    percy = prctile(x, [25 75]);
    hanni(3) = plotvertline(percy(1), [], '--', 'linewidth', 1, 'color', farbe);
    hanni(4) = plotvertline(percy(2), [], '--', 'linewidth', 1, 'color', farbe);
end
title(label)
    
if nargout > 0
    h = hanni;
end