function out = dcolred(dcolstr)
% BOLDCOLSTAR ... 
%  
%   ... 

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 09-Sep-2020 20:25:04 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 9.8.0.1451342 (R2020a) Update 5 
% FILENAME  : boldcolstar.m 



%% FUNCTION BODY [ boldcolstar.m ] 

% [pre, post] = strtok(dcolstr, '.');
% out = sprintf('{\\bf%s}.{\\bf%s}', pre, post(2:end));

out = strcat('\dcolcolor{darkred} ', dcolstr);
