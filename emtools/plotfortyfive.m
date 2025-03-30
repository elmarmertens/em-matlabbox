function h = plotfortyfive(varargin)
% PLOTFORTYFIVE ...
%
%   ...

%% VERSION INFO
% AUTHOR    : Elmar Mertens
% $DATE     : 18-Jul-2021 22:34:15 $
% $Revision : 1.00 $
% DEVELOPED : 9.10.0.1684407 (R2021a) Update 3
% FILENAME  : plotfortyfive.m

if isempty(varargin)
    h = plot(xlim, xlim, 'k--');
else
    h = plot(xlim, xlim, varargin{:});
end
