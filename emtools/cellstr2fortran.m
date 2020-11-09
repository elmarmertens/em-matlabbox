function cellstr2fortran(filename, x)
% CELLSTR2FORTRAN writes a cell of string values into a file row-by-row
%  
% Usage: cellstr2fortran(filename, x)
% 
% See also mat2fortran, vec2fortran

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 24-Jun-2012 17:25:16 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 7.7.0.471 (R2008b) 
% FILENAME  : vec2fortran.m 

fid = fopen(filename, 'wt');
fprintf(fid, '%s\n', x{:});
fclose(fid);