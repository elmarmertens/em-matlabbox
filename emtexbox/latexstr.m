function l = latexstr(s)
% function l = latexstr(s)
% converts a string into being "latex compatible"
% l = strrep(s, '_', '\_');

%   Coded by  Elmar Mertens, em@elmarmertens.com

l = s;
l = strrep(l, '_', '\_');
l = strrep(l, '&', '\&');
l = strrep(l, '%', '\%');

