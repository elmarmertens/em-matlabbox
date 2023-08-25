function [betadraw, resid, E, V] = bayesAR1SURdraw(y, ylag, Vresid, beta0, betaV0i, Ndraws, rndStream)
% bayesAR1SURdraw

%   Coded by  Elmar Mertens, em@elmarmertens.com

%% parse inputs
if nargin < 6 || isempty(Ndraws)
   Ndraws = 1;
end
if nargin < 7 || isempty(rndStream)
   rndStream = getDefaultStream;
end


%% draw SUR coefficients
Ny         = size(y,2);
Iy         = eye(Ny);

zdraws     = randn(rndStream, Ny, Ndraws);

H          = Iy / Vresid; 
ytilde     = y * H;

XX         = ylag' * ylag;
Vi         = betaV0i + XX .* H;
cholVi     = chol(Vi, 'upper');
cholinvVE  = transpose(cholVi) \ (betaV0i * beta0 + sum(ylag .* ytilde)');
betadraw   = cholVi \ (cholinvVE + zdraws);

%% additional output arguments
if nargout > 1 && (Ndraws == 1)
   resid       = y - ylag * betadraw';
else
   resid = [];
end

if nargout > 2
   E = cholVi \ cholinvVE;
end
if nargout > 3
   invcholV = cholVi \ Iy;
   V        = invcholV * transpose(invcholV);
end
