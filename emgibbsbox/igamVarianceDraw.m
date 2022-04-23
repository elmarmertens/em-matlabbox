function draw = igamVarianceDraw(resid, ssr0, dof0, dim)
% IGVARIANCEDRAW draw = igVarianceDraw(resid, ssr0, dof0);
%  
%   ... 

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 20-Apr-2022 12:36:23 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 9.12.0.1884302 (R2022a) 
% FILENAME  : igVarianceDraw.m 


if nargin < 4 || isempty(dim)
    dim = 1;
end

if isvector(resid)
    resid = resid(:);
    dim   = 1;
end

T            = size(resid, dim);
ssr          = sum(resid.^2, dim) + ssr0;
dof          = T + dof0;
% igamrnd      = @(alpha,beta) 1 ./ gamrnd(alpha, 1 ./ beta);
% draw         = igamrnd(dof * .5, ssr * .5);
draw         = 1 ./ gamrnd(dof * .5, 2 ./ ssr); 