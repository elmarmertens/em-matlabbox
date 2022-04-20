function draw = igamVarianceDraw(resid, ssr0, dof0)
% IGVARIANCEDRAW draw = igVarianceDraw(resid, ssr0, dof0);
%  
%   ... 

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 20-Apr-2022 12:36:23 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 9.12.0.1884302 (R2022a) 
% FILENAME  : igVarianceDraw.m 



%% FUNCTION BODY [ igVarianceDraw.m ]
ssr          = sum(resid(:).^2) + ssr0;
dof          = length(resid(:)) + dof0;
% igamrnd      = @(alpha,beta) 1 ./ gamrnd(alpha, 1 ./ beta);
% draw         = igamrnd(dof * .5, ssr * .5);
draw         = 1 ./ gamrnd(dof * .5, 2 ./ ssr); % use of ./ appears faster than /
