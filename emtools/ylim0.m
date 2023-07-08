function ylim0
% YLIM0 makes sure that that y-axcis includes 0
%
%   ...

%% VERSION INFO
% AUTHOR    : Elmar Mertens
% $DATE     : 10-May-2023 14:14:13 $
% $Revision : 1.00 $
% DEVELOPED : 9.14.0.2239454 (R2023a) Update 1
% FILENAME  : ylim0.m

YLIM = ylim;
YLIM(1) = min(0, YLIM(1));
YLIM(2) = max(0, YLIM(2));
ylim(YLIM)