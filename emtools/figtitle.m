function figtitle(s,fh)
% FIGTITLE ... 
%  
%   ... 

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 22-Feb-2021 23:12:13 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 9.9.0.1592791 (R2020b) Update 5 
% FILENAME  : figtitle.m 



if nargin < 2
    fh = gcf;
end

title(s)
set(fh, 'name', s)