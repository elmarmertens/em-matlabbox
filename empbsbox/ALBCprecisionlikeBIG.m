function [llf, CC, QQ, RR1, arows, acols, asortndx, brows, bcols, bsortndx] = ...
    ALBCprecisionlikeBIG(aaa,invbbb,ccc,y,x0,invsqrtsig0,CC,QQ,RR1,arows, acols, asortndx, brows, bcols, bsortndx)
% ALAGBCPRECISIONSAMPLER ...
%
% allows for lags of A; important: aaa should be ordered from p to 1 in 3rd dimension
%   ...

%% VERSION INFO
% AUTHOR    : Elmar Mertens

% get dimensions
[Ny, T] = size(y);
p       = size(aaa,3);
Nx      = size(aaa,1);
Nw      = size(invbbb,2);
if Nx ~= Nw
    error('dimension mismatch: Nx not equal to Nw')
end

if nargin < 7
    CC  = [];
    QQ  = [];
    RR1 = [];
    [arows, acols, asortndx, brows, bcols, bsortndx] = deal([]);
end

if ndims(aaa) <= 3
    aaa = repmat(aaa, [1 1 1 T]);
end
if ismatrix(invbbb)
    invbbb = repmat(invbbb, [1 1 T]);
end
if ismatrix(ccc)
    ccc = repmat(ccc, [1 1 T]);
end

Nx0   = Nx * p;
NyT   = Ny * T;
NxT   = Nx * T;
NxTp  = Nx * (T + p);


%% construct vectorized state space
Y     = reshape(y, NyT, 1);
XX0   = sparse(1:Nx0, 1, x0, NxTp, 1);

%% vectorize input matrices
NxNx         = Nx * Nx;
NxNxT        = NxNx * T;
invsqrtsig0  = reshape(invsqrtsig0, Nx0 * Nx0, 1); 
invbbb       = reshape(invbbb, NxNxT, 1);
ccc          = reshape(ccc, Ny * NxT, 1);

%% CC and prepare Arows and Brows

if isempty(CC)

    % AA
    arows1     = transpose(1 : NxTp);
    acols1     = transpose(1 : NxTp);

    arows2     = repmat((1 : Nx)', 1, Nx * p);
    arows2     = Nx0 + arows2 + permute(Nx * (0 : T - 1), [1 3 2]);
    acols2     = repmat(1 : Nx * p, Nx,1) + permute(Nx * (0 : T - 1), [1 3 2]);

    arows      = [arows1; reshape(arows2, NxNx * p * T, 1)];
    acols      = [acols1; reshape(acols2, NxNx * p * T, 1)];

    [acols, asortndx] = sort(acols);
    arows             = arows(asortndx);

    % BB
    brows0  = repmat((1 : Nx0)', 1 , Nx0);
    brows1  = Nx0 + repmat((1 : Nx)', 1 , Nx) + permute(Nx * (0 : T-1), [1 3 2]);
    brows   = [reshape(brows0, Nx0 * Nx0, 1); reshape(brows1, NxNx * T, 1)];

    bcols0  = repmat((1 : Nx0), Nx0, 1);
    bcols1  = Nx0 + repmat((1 : Nx), Nx, 1) + permute(Nx * (0 : T-1), [1 3 2]);
    bcols   = [reshape(bcols0, Nx0 * Nx0, 1); reshape(bcols1, NxNx * T, 1)];

    [bcols, bsortndx] = sort(bcols);
    brows             = brows(bsortndx);

    % C
    crows     = repmat((1 : Ny)', 1 , Nx, T) + permute(Ny * (0 : T-1), [1 3 2]);
    ccols     = Nx0 + repmat(1 : NxT, Ny, 1);
    CC        = sparse(reshape(crows, Ny * Nx * T, 1), reshape(ccols, Ny * Nx * T, 1), ccc, NyT, NxTp);

    % perform QR
    [QQ,RR]   = qr(CC');
    [N1, NN]  = size(CC);
    %     N2        = NN - N1;
    RR1       = RR(1:N1,1:N1)';

else

    N1        = size(RR1,1);
    NN        = size(QQ,1);
    %     N2        = NN - N1;
    
end

% QQ1       = QQ(:,1:N1)';
% QQ2       = QQ(:,N1+1:end)';

%% sparse builds for BB and AA
values  = [invsqrtsig0; invbbb];
values  = values(bsortndx);
invBB   = sparse(brows, bcols, values, NxTp, NxTp);

values1    = ones(NxTp,1);
values2    = reshape(-aaa, NxNx * p * T, 1); % (:,:,p:-1:1,:);
values     = [values1; values2];
values     = values(asortndx);
AA         = sparse(arows, acols, values, NxTp, NxTp);

%% means and innovations
EX        = AA \ XX0;
EY        = CC * EX;
% Ydev      = Y - EY;

% X1tilde   = RR1 \ (Y - EY);

%% likelihood calculations
AAtilde     = invBB * AA;
AAtildeQQ   = AAtilde * QQ;

Ptilde      = AAtildeQQ' * AAtildeQQ;

x1tilde     = RR1 \ (Y - EY);
logdetR11   = 2 * sum(log(abs(diag(RR1)))); % abs since RR1 is output from QR

%% tedious
sqrtPtilde     = chol(Ptilde, 'lower');
invsqrtPtilde  = eye(NN) / sqrtPtilde;
varXtilde      = invsqrtPtilde' * invsqrtPtilde;
sqrtvarX1tilde = chol(varXtilde(1:N1,1:N1), 'lower');
x1dev          = sqrtvarX1tilde \ x1tilde;
logdetX11      = 2 * sum(log(diag(sqrtvarX1tilde)));
llf            = -.5 * (N1 * log(2 * pi) + logdetR11 + logdetX11 + sum(x1dev.^2)); 


