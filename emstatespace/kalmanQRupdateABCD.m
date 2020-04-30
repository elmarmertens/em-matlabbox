function [sqrtOmega, K, sqrtSigma] = kalmanQRupdateABCD(A, B, C, D, sqrtSigmaLag)

[Nx,Nw] = size(B);
Ny      = size(C,1);

M = zeros(Ny + Nx, Nx + Nw);

% if Nx + Nw < Ny + Nx
%     error('QR requires rectangular matrix with more rows than cols') 
%     % recall that QR is applied to M'
% end

M(1:Ny, 1:Nx)           = C * sqrtSigmaLag;
M(1:Ny, Nx+(1:Nw))      = D;
M(Ny+(1:Nx), 1:Nx)      = A * sqrtSigmaLag;
M(Ny+(1:Nx), Nx+(1:Nw)) = B;


[~,R] = qr(M');

R = R';

sqrtOmega   = R(1:Ny,1:Ny);
K           = R(Ny+(1:Nx),1:Ny) / sqrtOmega;
sqrtSigma   = R(Ny+(1:Nx),Ny+(1:Nx));




