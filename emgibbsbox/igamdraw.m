function draws = igamdraw(ssr, dof, Ndraw)
% igamrnd 
%
% USAGE: igams = igamdraw(ssr, dof, Ndraw)
%
% See also: iwishdraw, igammadraw, gamrnd

%   Coded by  Elmar Mertens, em@elmarmertens.com

if ~isscalar(ssr)
    error('ssr is supposed to be scalar')
end
if ~isscalar(dof)
    error('dof is supposed to be scalar')
end
% igamrnd      = @(alpha,beta) 1 ./ gamrnd(alpha, 1 ./ beta);
draws  = 1 ./ gamrnd(dof * .5, 2 ./ ssr, Ndraw);
