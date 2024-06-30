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
   ndx = 1 : length(clrs);
else
   ndx(ndx>length(clrs)) = mod(ndx(ndx>length(clrs)), length(clrs));
   ndx(ndx == 0)         = length(clrs);
end

if isscalar(ndx)
    c = clrs{ndx};
else
    c = clrs(ndx);
end

