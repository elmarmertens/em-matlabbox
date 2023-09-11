function [YbarDraw, YgapDraw, Abar, Cbar, agaprows, agapcols, agapsortndx, bgaprows, bgapcols] = ...
    commontrendcyclePrecisisionsampler(y, agap, invbgap, invbbar, ybar0, rndStream, Abar, Cbar, ...
    agaprows, agapcols, agapsortndx, bgaprows, bgapcols)

%% VERSION INFO
% AUTHOR    : Elmar Mertens

if nargin < 7
    [Abar, agaprows, agapcols, agapsortndx, bgaprows, bgapcols] = deal([]);
end
% get dimensions
[Ny, T] = size(y);
p       = size(agap,3);

if ndims(agap) <= 3
    agap = repmat(agap, [1 1 1 T]);
end
if ismatrix(invbgap)
    invbgap = repmat(invbgap, [1 1 T]);
end
if ismatrix(invbbar)
    invbbar = repmat(invbbar, [1 1 T]);
end

NyT = Ny * T;

%% construct vectorized state space
Y      = reshape(y, NyT, 1);
Ybar0  = sparse(1, 1, ybar0, T, 1);

NyNy   = Ny * Ny;
Nagap  = NyT + NyNy * p * (T - p) + sum(NyNy * (1 : p - 1));
Nb     = NyNy * T;
    
if isempty(Abar)
    % Abar
    Abar = spdiags([-1 * ones(T,1), ones(T,1)], [-1, 0], T, T);

    % Agap
    [agaprows, agapcols] = deal(NaN(Nagap, 1));
    agaprows(1:NyT) = 1:NyT;
    agapcols(1:NyT) = 1:NyT;
    offset = NyT;
    for k = 1 : p
        theserows   = repmat((1 : Ny)', 1 , Ny, T - k);
        theserows   = theserows + permute(Ny * (k : T-1), [1 3 2]);

        thesecols  = repmat(1 : Ny * (T - k), Ny, 1);

        agaprows(offset + (1 : NyNy * (T-k))) = theserows(:);
        agapcols(offset + (1 : NyNy * (T-k))) = thesecols(:);

        offset = offset + NyNy * (T-k);
    end
    % sort Agap indices
    ndx = sub2ind([NyT, NyT], agaprows, agapcols);
    [~, agapsortndx] = sort(ndx);
    agaprows         = agaprows(agapsortndx);
    agapcols         = agapcols(agapsortndx);

    % Bgap 
    bgaprows = repmat((1:Ny)', 1, Ny) + permute(Ny * (0 : T - 1), [1 3 2]);
    bgapcols = repmat(1:Ny, Ny, 1) + permute(Ny * (0 : T - 1), [1 3 2]);
    bgaprows = reshape(bgaprows, Nb, 1);
    bgapcols = reshape(bgapcols, Nb, 1);
    % no sorting needed
    % ndx   = sub2ind([NyT, NyT], brows, bcols);
    % [~, bsortndx] = sort(ndx);
    % brows = brows(bsortndx);
    % bcols = bcols(bsortndx);

    % Cbar
    Ncbar = NyT;
    crows = 1:NyT;
    ccols = repmat(1:T, Ny, 1);
    Cbar  = sparse(crows(:), ccols(:), ones(Ncbar,1), NyT,T);
end

avalues           = ones(Nagap,1);
offset = NyT;
for k = 1 : p
    avalues(offset + (1 : NyNy * (T-k))) = -reshape(agap(:,:,k,1+k:T), NyNy * (T - k), 1);
    offset = offset + NyNy * (T-k);
end
avalues             = avalues(agapsortndx);
Agap                = sparse(agaprows, agapcols, avalues, NyT, NyT);

invBbar = spdiags(invbbar(:), 0, T, T);

invBgap = sparse(bgaprows, bgapcols, invbgap(:), NyT, NyT);

Abar      = invBbar * Abar;
P0bar     = Abar' * Abar;
Ybar0     = invBbar * Ybar0; % Ytilde0 is assumed zero

Agap      = invBgap * Agap * Cbar;
P0tilde   = Agap' * Agap;
Y         = invBgap * Y;

% posterior
Pbar         = P0bar + P0tilde;
sqrtPbar     = chol(Pbar, 'lower');
Pbarmu       = Abar' * Ybar0 + Agap' * Y;
sqrtPmu      = sqrtPbar \ Pbarmu;

Zdraw        = randn(rndStream, T, 1);
YbarDraw     = transpose(sqrtPbar) \ (sqrtPmu + Zdraw);
if nargout > 1
    YgapDraw     = Y - Cbar * YbarDraw;
end