function md = xls2mdate(xd, offset)
% XLS2MDATE ... 
%  
%   ... 

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 10-Feb-2021 18:07:31 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 9.9.0.1570001 (R2020b) Update 4 
% FILENAME  : xls2mdate.m 

if nargin < 2 
    offset = 693960; % 1900 system
    % offset = 695422; % 1904 system
end

md = xd + offset;
