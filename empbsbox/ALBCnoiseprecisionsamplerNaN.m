function [Xdraw,  XshockDraw, NoiseDraw, ...
    arows, acols, asortndx, brows, bcols, crows, ccols] = ...
    ALBCnoiseprecisionsamplerNaN(aaa,invbbb,ccc,invnoisevol,y,yNaN,x0,invsqrtsig0,rndStream, ...
    arows,acols,asortndx,brows,bcols,crows,ccols)

% ALBCnoiseprecisionsamplerNaN ...
%
% allows for lags of A; important: aaa should be ordered from p to 1 in 3rd dimension
%   ...
% note: invnoisevol should be Ny x T matrix or Ny * T vector, sampler assumes noise-vol is diagonal

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
if ismatrix(invnoisevol) && (size(invnoisevol,1) ~= Ny || size(invnoisevol,2) ~= T)
        error('dimension mismatch: invnoisevol should be Ny x T matrix (or Ny * T vector)')
elseif isvector(invnoisevol) && length(invnoisevol) ~= Ny * T
    error('dimension mismatch: invnoisevol should be Ny x T matrix (or Ny * T vector)')
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
Ynan  = reshape(yNaN, NyT, 1);
Y     = Y(~Ynan);
if length(x0) == NxTp
    XX0 = x0;
else % for now assume x0 has Nx0=Nx*p length unless it has full length
    XX0 = sparse(1:Nx0, 1, x0, NxTp, 1); 
end

%% vectorize input matrices
NxNx         = Nx * Nx;
NxNxT        = NxNx * T;
invsqrtsig0  = reshape(invsqrtsig0, Nx0 * Nx0, 1);
invbbb       = reshape(invbbb, NxNxT, 1);
ccc          = reshape(ccc, Ny * NxT, 1);

%% construct row- and column indices (if needed)
if nargin < 10
    % AA
    arows1     = transpose(1 : NxTp);
    acols1     = transpose(1 : NxTp);

    arows2     = repmat((1 : Nx)', 1, Nx * p);
    arows2     = Nx0 + arows2 + permute(Nx * (0 : T - 1), [1 3 2]);
    acols2     = repmat(1 : Nx * p, Nx,1) + permute(Nx * (0 : T - 1), [1 3 2]);

    arows      = [arows1; reshape(arows2, NxNx * p * T, 1)];
    acols      = [acols1; reshape(acols2, NxNx * p * T, 1)];

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

    % sort B indices -- not necessary, since the above is already sorted
    % ndx = sub2ind([NxTp, NxTp], brows, bcols);
    % [~, bsortndx] = sort(ndx);
    % brows         = brows(bsortndx);
    % bcols         = bcols(bsortndx);

    %% CC
    crows     = repmat((1 : Ny)', 1 , Nx, T) + permute(Ny * (0 : T-1), [1 3 2]);
    ccols     = Nx0 + repmat(1 : NxT, Ny, 1);
    crows     = crows(:);
    ccols     = ccols(:);


    % sort C indices -- not necessary, since the above is already sorted
    % ndx = sub2ind([NyT, NxTp], crows, ccols);
    % if ~issorted(ndx)
    %     error houston
    % end
    % [~, csortndx] = sort(ndx);
    % crows         = crows(csortndx);
    % ccols         = ccols(csortndx);

end
%% CC and prepare Arows and Brows

% AA
% values1    = ones(NxTp,1);
% values2    = reshape(-aaa, NxNx * p * T, 1); % (:,:,p:-1:1,:);
% values     = [values1; values2];
values             = ones(NxTp + NxNx * p * T,1);
values(NxTp+1:end) = -aaa(:);
values             = values(asortndx);
AA                 = sparse(arows, acols, values, NxTp, NxTp);

% BB
values        = [invsqrtsig0; invbbb];
invsqrtSIGMA  = sparse(brows, bcols, values, NxTp, NxTp);

% C
CC = sparse(crows, ccols, ccc, NyT, NxTp);
% drop rows associated with NaN
CC = CC(~yNaN,:);
	
% noiseVOL
invsqrtOMEGA = spdiags(invnoisevol(:),0,NyT,NyT);
invsqrtOMEGA = invsqrtOMEGA(~yNaN,~yNaN);


%% means and innovations
AAtilde            = invsqrtSIGMA * AA;
XX0tilde           = invsqrtSIGMA * XX0;

CCtilde            = invsqrtOMEGA * CC;
Ytilde             = invsqrtOMEGA * Y;

P                   = AAtilde' * AAtilde + (CCtilde' * CCtilde);
[sqrtP, flag]       = chol(P, 'lower');

if flag > 0
    warning('P not posdf, using QR instead')
    % via qr -- much slower
    M = [AAtilde; CCtilde];
    m = size(M,2);
    [~, R] = qr(M);
    sqrtP = R(1:m,1:m)';

    % checkdiff(sqrtP * sqrtP', P);
    % checkdiff(M'*M,P);
end

sqrtPXhat    = sqrtP \ (AAtilde' * XX0tilde + CCtilde' * Ytilde);
Zdraw        = randn(rndStream, NxTp, 1);
Xdraw        = transpose(sqrtP) \ (sqrtPXhat + Zdraw);
if nargout > 1
    XshockDraw   = AA * Xdraw - XX0;
end
if nargout > 2
    NoiseDraw    = Y - CC * Xdraw;
end

