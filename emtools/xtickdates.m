function xtickdates(dates, varargin)
% XTICKDATES sets dateticks 
%  xtickdates(dates) or xtickdates(dates, ...) set date limits to conform to date vector dates
%  and passes any additional arguments to Matlab's datetick function
%
% see also datetick

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 06-Dec-2012 10:38:58 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 7.14.0.739 (R2012a) 
% FILENAME  : xtickdates.m 

xlim(dates([1 end]))
datetick('x', 'keeplimits', varargin{:})
% set(gcf, 'Renderer', 'painters') % legacy
