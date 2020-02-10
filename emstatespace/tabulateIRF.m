function tabulateIRF(irf, ynames, wnames)
% function tabulateIRF(irf, ynames, wnames)
% print out each IRF on the screen

[ny, nw, lags] = size(irf);
lags = lags - 1;

for y = 1 : ny
   hrulefill
   fprintf('IRF of %s\n', ynames{y});
   fprintf('%10s \t', 'lags', wnames{:});
   fprintf('\n');
   for t = 0 : lags
      fprintf('%10.4f \t', t, squeeze(irf(y,:,t+1)))
      fprintf('\n');
   end
end
fprintf('\n');