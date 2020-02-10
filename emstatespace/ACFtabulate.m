function TAB = ACFtabulate(ACF, k, n)
% function TAB = ACFtabulate(ACF, k, n)
% Tabulates ACF (+/- k lags) of every element with the n'th element
% note: this table is in the spirit of the HP tables, e.g. in the Cooley/Prescott book

% Elmar Mertens
% www.elmarmertens.ch

error(nargchk(3,3,nargin));

[N, N1, K, jack] = size(ACF);

if jack > 1
   error('ACF should be 3D')
end

if isempty(k)
    k = K - 1;
end
% if nargin < 3
%     k = K - 1;
%     n = k;
% end

if N ~= N1
    error('dimension mismatch')
end

TAB = repmat(NaN, [N, 2 * k + 1]);

TAB(:, (k + 1)) = ACF(:, n, 1);

for i = 1 : k % loop over leads and lags
    TAB(:, (k + 1) + i) = ACF(n, :, i + 1)';
    TAB(:, (k + 1) - i) = ACF(:, n, i + 1);
end