function ACF = data2ACF(Y, K, SameSampleFlag, dofscale, demeanflag)
% function ACF = data2ACF(Y, K, SameSampleFlag, dofscale, demeanflag)
% defaults:
% K               = size(Y,2)-1
% SameSampleFlag  = false
% dofscale        = false
% demeanflag      = true
% NOTE: default setting yield sample ACF

% note: acf2bartlettspectrm0(data2ACF(Y,k,false,false), k) reproduces haccme(Y',k)

error(nargchk(1,5,nargin))

[T, N] = size(Y);

if nargin < 2 || isempty(K)
    K = T-1;
end
if nargin < 3 || isempty(SameSampleFlag)
    SameSampleFlag = false;
end
if nargin < 4 || isempty(dofscale)
   dofscale = false;
elseif dofscale && SameSampleFlag
   warning('elmi:msg', 'SameSampleFlag is set. Ignoring dofscale == true')
end
if nargin < 5 || isempty(demeanflag)
    demeanflag = true;
end


if SameSampleFlag
    Y = Y(K + 1 : end, :);
end

if demeanflag
    Y = bsxfun(@minus, Y, mean(Y));
end

% ACF = repmat(NaN, [N, N, K + 1]);
ACF = zeros(N, N, K + 1);

for k = 0 : K
   if SameSampleFlag
      
      ACF(:, :, k + 1) = Y(K + 1 : end, :)' * Y((K + 1: end) - k, :) / (T-K); 
      
   else
      
      ACF(:, :, k + 1) = Y(k + 1 : end, :)' * Y(1 : end - k, :);
      
      if dofscale
         ACF(:, :, k + 1) = ACF(:, :, k + 1) / (T -k);
      else
         ACF(:, :, k + 1) = ACF(:, :, k + 1) / T;
      end
      
   end
end

