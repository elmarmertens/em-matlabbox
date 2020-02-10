function ACC = ACF2ACC(ACF, lowmemory)
% function ACC = ACF2ACC(ACF, lowmemory)
% converts autocovariance function ACF into autocorrelation function ACC
% first three dimensions must be: nvar x nvar x lags
% works with arbitrary dimensions beyond 3 (e.g.: monte carlo of conditional acf would be stored in a 5D ACF)
%
% "lowmemory = true" is more memory efficient, but uses more loops
% default: lowmemory = false

%   Coded by  Elmar Mertens, em@elmarmertens.com

% Elmar Mertens
% www.elmarmertens.ch

if nargin < 2
   lowmemory = false;
end


N = size(ACF, 1);
if N ~= size(ACF, 2);
   error('First two dimensions of ACF should be identical! size(ACF. 1:2) = [%d %d]', N, size(ACF, 2))
end

acf4D = reshape(ACF, N, N, size(ACF, 3), []); % collapse all dimensions above 4 into 4th dimension
ACC   = zeros(size(acf4D));
K     = size(acf4D, 3);
four  = size(acf4D, 4);
if lowmemory
   oo = zeros(N);
else
   oo = zeros(N, N, size(acf4D, 3));
end

for f = 1 : four

   if lowmemory

      vols = sqrt(diag(acf4D(:,:,1,f)));
      vv   = (vols * vols');
      zndx = vv ~= 0;

      for k = 1 : K
         jim         = acf4D(:, :, k,f);
         jack        = oo;
         jack(zndx)  = jim(zndx) ./ vv(zndx);
         jack(~zndx) = Inf;

         ACC(:, :, k,f) = jack;
      end


   else % good memory conditions, replace the k-loop with a repmat

      vols              = sqrt(diag(acf4D(:,:,1,f)));
      vols(vols == 0)   = NaN; % to avoid divide-by-zero warnings
      jack              = bsxfun(@rdivide,acf4D(:,:,:,f),vols);
      ACC(:, :, :,f)    = bsxfun(@rdivide,jack,vols');

   end


end

% convert back to original dimensions
ACC = reshape(ACC, size(ACF));


%% move rounding errors inside unit circle and issue warning about larger errors
outside = abs(ACC) > 1;
if all((abs(ACC(outside)) - 1) < 1e-14)
   ACC(outside) = sign(ACC(outside));
else
   warning('em:base', 'Some ACC significantly outside unit circle, max deviation: %12.6e', max(abs(ACC(outside)) - 1))
end
