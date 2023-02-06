function writedatatable2tex(wrap, tabname, dates, data, datalabels, datefmt)
% WRITEDATATABLE stores times series data in a csv file with column headers
% and a date column.
%
% USAGE: writedatatable(wrap, filename, dates, data, datalabels, datefmt)
% - wrap determines the directory in whic the csv file is to be stored. wrap
%   can be one of the following:
%   -- empty: the current directory (pwd) will be used)
%   -- a string, designating the path to the target directory
%   - a structure with a field called "dir" designating the path to the target directory.
% - filename: name of the target file (without "csv" suffix, which will be appended)
% - dates: Tx1 vector of matlab date numbers
% - data: TxN data matrix (can contain NaNs)
% - datalabels: cell of length N, containing colun labels 
% - datefmt: matlab formatting string for outputting the dates in column 1 
%   (default: % 'yyyyqq')
%  
%   See also: dlmwrite, csvwrite, pwd, datestr

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 28-Mar-2016 10:20:04 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 8.6.0.267246 (R2015b) 
% FILENAME  : writedatatable.m 

if nargin < 6
    datefmt = 'yyyyqq';
end


if isempty(wrap) 
    return
else
    if isstruct(wrap) && isfield(wrap, 'dir')
        tabdir = wrap.dir;
    else % assuming it is a string ....
        tabdir = wrap;
    end
end

fid = fopen(fullfile(tabdir, strcat(tabname, '.tex')), 'wt');

fprintf(fid, '\\begin{small}\n');
fprintf(fid, '\\begin{center}\n');
fprintf(fid, '\\begin{tabular}{l%s}\n', repmat('.4', 1, size(data,2)));
fprintf(fid, '\\toprule\n');
fprintf(fid, 'Dates ');
fprintf(fid, ' & \\multicolumn{1}{c}{%s}', datalabels{:});
fprintf(fid, '\\\\\n');
fprintf(fid, '\\midrule\n');
for t = 1 : size(data,1)
     fprintf(fid, '%s ', datestr(dates(t), datefmt)); %#ok<DATST> 
     fprintf(fid, '& %6.2f ', data(t,:));
     fprintf(fid, '\\\\\n');
end
 
fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{center}\n');
fprintf(fid, '\n');
fprintf(fid, '\\end{small}\n');
fclose(fid);
type(fullfile(tabdir, strcat(tabname, '.tex')))
latexwrapper(wrap, 'add', 'sidetab', strcat(tabname, '.tex'), [])
