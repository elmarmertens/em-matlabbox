function [PAI, sqrtXX, sqrtSSR, dof, R2, R2yx, R2yz, R2zx] = BVARZjeffries(Y,X,Z)
% BVARZJEFF BVAR estimation with jeffries prior and controls Z, orthogonalized parametrization
%
% Y = X * PAI + Ztilde * GAMMA + E
% where Ztilde is Z orthogonalized against X (residual of regressing Z on X)
%
% Note: PAI is the coefficient on X in the Ztilde parametrization, which equals
% the coefficient from regressing Y on X alone; it is NOT the coefficient on X
% in a joint regression of Y on X and raw Z. sqrtSSR is the SSR of the joint
% regression, and dof = T - (Nx + Nz).
%
% See also: BVARjeffries, BVARdraws

% get dimensions
[T, Nx]  = size(X);
[~, Nz]  = size(Z);
[~, Ny]  = size(Y);

K = Nx + Nz;

if T < K + Ny
    error('em:BVARZjeffries', 'need T >= Nx + Nz + Ny (T=%d, Nx=%d, Nz=%d, Ny=%d)', T, Nx, Nz, Ny)
end
% perform OLS via QR decomposition
RR      = qr([X Z Y],'econ'); % single output skips forming Q
RR      = RR';

Xndx = 1:Nx;
Zndx = Nx+(1:Nz);
Yndx = Nx+Nz+(1:Ny);
% select output
sqrtXX    = RR(Xndx,Xndx);
sqrtZX    = RR(Zndx,Xndx);
% sqrtZZ    = RR(Zndx,Zndx);

sqrtYX    = RR(Yndx,Xndx);
sqrtYZ    = RR(Yndx,Zndx);

sqrtSSR   = RR(Yndx,Yndx);

PAI   = (sqrtYX / sqrtXX)';
% GAMMA = (sqrtYZ / sqrtZZ)';

dof = T - K;

% report R2
if nargout > 4
    YY   = sum(Y.^2, 1)';
    R2   = 1 - sum(sqrtSSR.^2, 2) ./ YY;
    R2yz = sum(sqrtYZ.^2, 2) ./ YY;
    R2yx = sum(sqrtYX.^2, 2) ./ YY;
    R2zx = sum(sqrtZX.^2, 2) ./ sum(Z.^2, 1)';
end

% checkdiff(sqrtXX * sqrtXX', X' * X);
% checkdiff(sqrtYX * sqrtYX' + sqrtYZ * sqrtYZ' + sqrtSSR * sqrtSSR', Y' * Y);
