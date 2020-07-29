function int2fortran(filename, x)
% INT2FORTRAN writes a matrix of integerns in a text file for import into fortran 
%  
% Usage: int2fortran(filename, x)
% 
% See also mat2fortran, logical2fortran

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 24-Jun-2012 17:25:16 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 7.7.0.471 (R2008b) 
% FILENAME  : vec2fortran.m 

dlmwrite(filename, x, 'delimiter','', 'precision', '%10d');