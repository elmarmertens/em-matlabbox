function drps = drpsDraws(y, Xdraws, binEdges,  weights, hasBeenSorted)
% DRPSDRAWS ...
%
% computes DRPS, allows to handle 2D inputs
%   ...

% TODO: hasBeenSorted not yet tested

%% check inputs
narginchk(3,5)

% X
if ~ismatrix(Xdraws) 
    error('ND Xdraws not supported')
end
[Ndraws, Nvar] = size(Xdraws);
% y
y = transpose(y(:)); % y is a row vector from now on
if Nvar ~= length(y)
    error('dimension mismatch Nvar')
end
% NaNs
yNaN = isnan(y);
if any(isnan(Xdraws), 'all') || all(yNaN, 'all')
    drps  = NaN(1,Nvar);
    return;
end

binEdges = binEdges(:);
Nbins    = length(binEdges);
binEdges = permute(binEdges, [2 3 1]); % 1 x 1 x Nbins from now on

if nargin < 4 || isempty(weights)
    weights = ones(Ndraws,1) ./ Ndraws;
else
    if any(abs(sum(weights) - 1) > 1e-6, 'all')
        error('weights do not sum to one')
    end
end

if nargin < 5 
    hasBeenSorted = [];
end
if isempty(hasBeenSorted)
    hasBeenSorted = issorted(Xdraws);
end


%% sort if needed
if hasBeenSorted
    Xordered = Xdraws;
    weights  = repmat(weights, [1 Nvar Nbins]);
else
    [Xordered, sortNdx] = sort(Xdraws, 1); % sort columnwise
    if ~isvector(weights)
        error('weights should be vector')
    end
    weights = weights(:); % ensure it is a column vector
    weights = weights(sortNdx); % weights corresponding to sorted draws, automatically expanded to Ndraws x Nvar matrix
    weights = repmat(weights, [1 1 Nbins]); % repeated for each of the bins
end

%% compute DRPS
Yhit            = permute(y < binEdges, [2 3 1]); 
Yhit            = double(Yhit);
Yhit(yNaN,:)    = NaN;

ndxCDF          = Xordered > binEdges; % 1 if a draw is greater than a given bin edge, o/w zero
weights(ndxCDF) = 0; % weights of draws corresponding to draws larger than a given bin edge are zeroed out
binCDF          = permute(sum(weights, 1), [2 3 1]); % using squeeze would problematic if Nvar = 1

drps            = sum((binCDF - Yhit).^2, 2);




