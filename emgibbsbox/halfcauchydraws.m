function draws = halfcauchydraws(varargin)
% HALFCAUCHYDRAWS ... 
% 
% draws = halfcauchydraws(Ndraws)
%   ... 

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 23-Apr-2022 16:31:33 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 9.12.0.1884302 (R2022a) 
% FILENAME  : halfcauchydraws.m 


draws      = abs(trnd(1, varargin{:}));
    

% oneN2      = ones(N,2);
% gdraws     = gamrnd(.5, oneN2);
% draws      = gdraws(:,1) ./ gdraws(:,2);

% faster; in particular for N below 1e4
