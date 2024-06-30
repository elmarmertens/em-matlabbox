function xlim0
% XLIM0 makes sure that that y-axis includes 0
%
%   ...

%% VERSION INFO
% AUTHOR    : Elmar Mertens
% $DATE     : 10-May-2023 14:14:13 $
% $Revision : 1.00 $
% DEVELOPED : 9.14.0.2239454 (R2023a) Update 1
% FILENAME  : xlim0.m

XLIM = xlim;
XLIM(1) = min(0, XLIM(1));
XLIM(2) = max(0, XLIM(2));
xlim(XLIM)