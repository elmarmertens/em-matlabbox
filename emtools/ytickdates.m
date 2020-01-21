function ytickdates(dates, varargin)
% YTICKDATES sets dateticks 
%  ytickdates(dates) or ytickdates(dates, 10, 'keepticks') etc
%  in addition to calling datetick, the function also implements a bugfix
%  to ensure properly prionted datelabels by disablign 3GL rendering
% see also datetick

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 06-Dec-2012 10:38:58 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 7.14.0.739 (R2012a) 
% FILENAME  : xtickdates.m 

ylim(dates([1 end]))
datetick('y', 'keeplimits', varargin{:})
% set(gcf, 'Renderer', 'painters')
