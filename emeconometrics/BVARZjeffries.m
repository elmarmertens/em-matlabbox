function [PAI, sqrtXX, sqrtSSR, dof, R2, R2yx, R2yz, R2zx] = BVARZjeffries(Y,X,Z)
% BVARZJEFF ... 
%  
% Y = X * PAI + Ztilde * GAMMA + E
%   ... 

% get dimensions
[T, Nx]  = size(X);
[~, Nz]  = size(Z);
[~, Ny]  = size(Y);

K = Nx + Nz;
% perform OLS via QR decomposition
[~,RR]  = qr([X Z Y],0);
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
    R2   = 1 - diag(sqrtSSR * sqrtSSR') ./ diag(Y' * Y);
    R2yz = diag(sqrtYZ * sqrtYZ') ./ diag(Y' * Y);
    R2yx = diag(sqrtYX * sqrtYX') ./ diag(Y' * Y);
    R2zx = diag(sqrtZX * sqrtZX') ./ diag(Z' * Z);
end

% checkdiff(sqrtXX * sqrtXX', X' * X);
% checkdiff(sqrtYX * sqrtYX' + sqrtYZ * sqrtYZ' + sqrtSSR * sqrtSSR', Y' * Y);
