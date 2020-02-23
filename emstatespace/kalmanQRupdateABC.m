function [sqrtOmega, Ktilde, sqrtSigma] = kalmanQRupdateABC(A, B, C, sqrtSigmaLag)

[Nx,Nw] = size(B);
Ny      = size(C,1);

M = zeros(Ny + Nx, Ny + Nx + Nw);

M(Ny+(1:Nx), Ny+(1:Nx))    = A * sqrtSigmaLag;
M(Ny+(1:Nx), Ny+Nx+(1:Nw)) = B;
M(1:Ny,Ny+1:end)          = C * M(Ny+1:end,Ny+1:end);

[~,R] = qr(M');

R = R';

sqrtOmega   = R(1:Ny,1:Ny);
Ktilde      = R(Ny+(1:Nx),1:Ny);
sqrtSigma   = R(Ny+(1:Nx),Ny+(1:Nx));




