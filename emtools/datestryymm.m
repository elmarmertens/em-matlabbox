function s = datestryymm(d)
% DATESTRYYMM returns date number as string in format yyyy:Mmm 
%  
% USAGE: s = datestryymm(d)
%   ... 

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 28-Jul-2020 19:59:11 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 9.8.0.1417392 (R2020a) Update 4 
% FILENAME  : datestryymm.m 



s = sprintf('%s:M%s', datestr(d, 'yyyy'), datestr(d, 'mm'));