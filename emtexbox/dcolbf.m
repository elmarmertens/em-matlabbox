function out = dcolbf(x, fmt, boldflag, framenum)
% BOLDCOLSTAR ...
%
%   ...

%% VERSION INFO
% AUTHOR    : Elmar Mertens
% $DATE     : 09-Sep-2020 20:25:04 $
% $Revision : 1.00 $
% DEVELOPED : 9.8.0.1451342 (R2020a) Update 5
% FILENAME  : boldcolstar.m

if nargin < 2 || isempty(fmt)
    fmt = '%6.4f';
end
if nargin < 3 || isempty(boldflag)
    dobold = false;
else
    if isa(boldflag,'function_handle')
        dobold = boldflag(x);
    else
        dobold = boldflag;
    end
end
if nargin < 4
    framenum = [];
end
fmtstr = sprintf(fmt, x);
if dobold
    [pre, post] = strtok(fmtstr, '.');
    if isempty(framenum)
        out = sprintf('{\\textbf{%s}}.{\\textbf{%s}}', pre, post(2:end));
    else
        out = sprintf('{\\textbf<%d>{%s}}.{\\textbf<%d>{%s}}', framenum, pre, framenum, post(2:end));
    end
else
    out = fmtstr;
end


