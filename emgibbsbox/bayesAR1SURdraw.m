function [rhodraw, resid, Erho, Vrho] = bayesAR1SURdraw(y, ylag, Vresid, rho0, rhoV0i, Ndraws, rndStream)
% bayesAR1SURdraw

%   Coded by  Elmar Mertens, em@elmarmertens.com


if nargin < 6 || isempty(Ndraws)
   Ndraws = 1;
end
if nargin < 7 || isempty(rndStream)
   rndStream = getDefaultStream;
end


[T, Ny] = size(y);
resid   = NaN(T,Ny);
% rhodraw = NaN(Ny,1);
Iy      = eye(Ny);

H      = Iy / Vresid;
ytilde = y * H;
XX     = ylag' * ylag;

rhoVi = rhoV0i + XX .* H;
Vrho  = Iy / rhoVi;
Erho  = Vrho * (rhoV0i * rho0 + sum(ylag .* ytilde)');

cholVrho = chol(Vrho)';
rhodraw = bsxfun(@plus, Erho, cholVrho * randn(rndStream, Ny, Ndraws)); 

if nargout > 1 && (Ndraws == 1)
   resid       = y - bsxfun(@times, ylag, rhodraw');
end
   
