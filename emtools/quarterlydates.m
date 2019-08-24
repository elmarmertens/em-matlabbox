function qdates = quarterlydates(dates)
% QUARTERLYDATES converts dates into quarters 
%
%   qdates = quarterlydates(mdates)
% 
% returns a vector qdates which is of the same size as dates
% the datenumbers qdates correspond to the first day of each quarter 
%
% See also genrQdates, genrMdates


%   Coded by  Elmar Mertens, em@elmarmertens.com

[Y, M]  = datevec(dates);
Q       = floor((M - 1) / 3) + 1;
MQ      = (Q - 1) * 3 + 1;
qdates  = datenum(Y, MQ, 1);


