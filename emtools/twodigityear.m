function Y = twodigityear(d)
% function Y = year(d);

%   Coded by  Elmar Mertens, em@elmarmertens.com

jack = datevec(d);
Y    = jack(:,1);

Y    = rem(Y, 100);