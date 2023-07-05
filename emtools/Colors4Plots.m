function c = Colors4Plots(ndx)
% function c = Colors4Plots(ndx)
% clrs = {
% [0      0.4470 0.7410]; % 1: dark blue
% [0.8500 0.3250 0.0980]; % 2: orange
% [0.9290 0.6940 0.1250]; % 3: dark yellow
% [0.4940 0.1840 0.5560]; % 4: dark purple
% [0.4660 0.6740 0.1880]; % 5: medium green
% [0.3010 0.7450 0.9330]; % 6: light blue
% [0.6350 0.0780 0.1840]  % 7: dark red
% [0      0      0     ]  % 8: black
% }
%
% c = clrs(ndx);

%   Coded by  Elmar Mertens, em@elmarmertens.com


persistent clrs

clrs = {
    [0 0.4470 0.7410];      % dark blue
    [0.8500 0.3250 0.0980]; % orange
    [0.9290 0.6940 0.1250]; % dark yellow
    [0.4940 0.1840 0.5560]; % dark purple
    [0.4660 0.6740 0.1880]; % medium green
    [0.3010 0.7450 0.9330]; % light blue
    [0.6350 0.0780 0.1840]; % dark red
    [0      0      0     ]  % black
    };


if nargin < 1
    c = clrs;
else
    if ~isnumeric(ndx)
        s = ndx;
        switch lower(s)
            case {'darkblue', 'dark blue', 'blue'}
                ndx = 1;
            case {'orange'}
                ndx = 2;
            case {'dark yellow', 'darkyellow', 'yellow'}
                ndx = 3;
            case {'dark purple', 'darkpurple', 'purple'}
                ndx = 4;
            case {'medium green', 'mediumgreen', 'green'}
                ndx = 5;
            case {'light blue', 'lightblue'}
                ndx = 6;
            case {'dark red', 'darkred', 'red'}
                ndx = 7;
            otherwise % black
                ndx = 8;
        end
    end
    if isscalar(ndx)
        c = clrs{ndx};
    else
        c = clrs(ndx);
    end
end