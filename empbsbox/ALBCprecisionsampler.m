function [Xdraw, XshockDraw, CC, QQ, RR1, ...
    arows, acols, a0ndx, asortndx, brows, bcols, b0ndx, bsortndx] = ...
    ALBCprecisionsampler(aaa,invbbb,ccc,y,x0,invsqrtsig0,rndStream,CC,QQ,RR1,...
    arows, acols, a0ndx, asortndx, brows, bcols, b0ndx, bsortndx)
	
% ALBCprecisionsampler ...
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

if nargin < 8
    CC  = [];
    QQ  = [];
    RR1 = [];
    [arows, acols, a0ndx, asortndx, brows, bcols, b0ndx, bsortndx] = deal([]);
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
% ccc          = reshape(ccc, Ny * NxT, 1); % will be vectorized later

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

    % a0ndx
    avalues                 = ones(NxTp + NxNx * p * T,1);
    avalues(NxTp + 1 :end)  = -aaa(:);
    a0ndx                  = ~(avalues == 0);
    arows                  = arows(a0ndx);
    acols                  = acols(a0ndx);

    % sort A indices
    ndx = sub2ind([NxTp, NxTp], arows, acols);
    [~, asortndx] = sort(ndx);
    arows         = arows(asortndx);
    acols         = acols(asortndx);

    % BB
    brows0  = repmat((1 : Nx0)', 1 , Nx0);
    brows1  = Nx0 + repmat((1 : Nx)', 1 , Nx) + permute(Nx * (0 : T-1), [1 3 2]);
    brows   = [reshape(brows0, Nx0 * Nx0, 1); reshape(brows1, NxNx * T, 1)];

    bcols0  = repmat((1 : Nx0), Nx0, 1);
    bcols1  = Nx0 + repmat((1 : Nx), Nx, 1) + permute(Nx * (0 : T-1), [1 3 2]);
    bcols   = [reshape(bcols0, Nx0 * Nx0, 1); reshape(bcols1, NxNx * T, 1)];

    % collect zero indices
    bvalues  = [invsqrtsig0; invbbb];
    b0ndx   = ~(bvalues == 0);
    brows   = brows(b0ndx);
    bcols   = bcols(b0ndx);

    % sort B indices
    ndx = sub2ind([NxTp, NxTp], brows, bcols);
    [~, bsortndx] = sort(ndx);
    brows         = brows(bsortndx);
    bcols         = bcols(bsortndx);

    %% CC
    crows     = repmat((1 : Ny)', 1 , Nx, T) + permute(Ny * (0 : T-1), [1 3 2]);
    ccols     = Nx0 + repmat(1 : NxT, Ny, 1);
    crows     = crows(:);
    ccols     = ccols(:);


    % sort C indices
    ndx = sub2ind([NyT, NxTp], crows, ccols);
    [~, csortndx] = sort(ndx);
    crows         = crows(csortndx);
    ccols         = ccols(csortndx);
    ccc           = ccc(csortndx);
    CC            = sparse(reshape(crows, Ny * Nx * T, 1), reshape(ccols, Ny * Nx * T, 1), ccc, NyT, NxTp);

    % perform QR
    [QQ,RR]   = qr(CC');
    [N1, N2]  = size(CC);
    N2        = N2 - N1;
    RR1       = RR(1:N1,1:N1)';

else

    N1        = size(RR1,1);
    N2        = size(QQ,1) - N1;

    avalues                 = ones(NxTp + NxNx * p * T,1);
    avalues(NxTp + 1 :end)  = -aaa(:);

    bvalues  = [invsqrtsig0; invbbb];

end

QQ1       = QQ(:,1:N1)';
QQ2       = QQ(:,N1+1:end)';

%% sparse builds for BB and AA
bvalues  = bvalues(b0ndx);
bvalues  = bvalues(bsortndx);
invBB    = sparse(brows, bcols, bvalues, NxTp, NxTp);

avalues  = avalues(a0ndx);
avalues  = avalues(asortndx);
AA       = sparse(arows, acols, avalues, NxTp, NxTp);

%% means and innovations
EX        = AA \ XX0;
EY        = CC * EX;

X1tilde   = RR1 \ (Y - EY);

QQX1tilde = QQ1' * X1tilde;

%% precision-based sampler
AAtilde       = invBB * AA;
AAtildeQQX1   = AAtilde * QQX1tilde;
AAtildeQQ2    = AAtilde * QQ2';
invQSIG22     = transpose(AAtildeQQ2) * AAtildeQQ2;
cholinvQSIG22 = chol(invQSIG22, 'lower');

X2hat         = - cholinvQSIG22 \ (AAtildeQQ2' * AAtildeQQX1);

Z2draw        = randn(rndStream, N2, 1) + X2hat;
X2draw        = cholinvQSIG22' \ Z2draw;
Xdraw         = EX + QQX1tilde + QQ2' * X2draw;
if nargout > 1
    XshockDraw   = AA * Xdraw - XX0;
end

