function h = plotACFC(ACF, k, VariableNames, ACFfrac, ACCfrac, tickSize, sameLimFlag)
% PLOTACFC
% plot ACC on lower and ACF on upper diagonal
% if ACFfrac and ACCfrac can be used to superimpose additional ACF/ACC onto ACF plot:
% - grey shaded confidence intervals if they are matrices (4D, with percentiles in 4th dimension)
% - simple lines if they are cells (containing extra ACF)
%
% ACFfrac, and ACCfrac need to be same dimension as ACF and can contain percentiles of simulated ACF/ACC distribution
%
% USAGE: plotACFC(ACF, k, VariableNames, ACFfrac, ACCfrac, tickSize, sameLimFlag)

error(nargchk(1, 7, nargin))

ACC = ACF2ACC(ACF);

[N, N1, K] = size(ACF);
if N ~= N1
   error('dimension mismatch')
end

if nargin < 2 || isempty(k)
   k = K - 1;
end

if nargin < 3
   VariableNames = cell(N,1);
end

if nargin < 4
   ACFfrac = [];
end
if iscell(ACFfrac)
   otherACF = ACFfrac;
   ACFfrac = [];
   for o = 1 : length(otherACF)
      if any(size(ACF) ~= size(otherACF{o}))
         error('Dimension mismatch amongst the multiple ACF''s')
      end
      otherACC{o} = ACF2ACC(otherACF{o});
   end
else
   otherACF = [];
   otherACC = [];
end

if nargin < 5
   ACCfrac = [];
end

if nargin < 6 || isempty(tickSize)
   tickSize = 4; % min(k, floor(40 / (2 * k) * N));
end

if nargin < 7
   sameLimFlag = false;
end

% finished checking inputs

% determine min/maxy (disregard other*)
jockey = ACF(:,:,1:k+1);
miny = min(jockey(:));
maxy = max(jockey(:));

if ~isempty(ACFfrac)
   jockey = ACFfrac(:, :, 1:k+1, :);
   miny = min(miny, min(jockey(:)));
   maxy = max(miny, max(jockey(:)));
end

% start plotting
h = zeros(N,N);
figure
for c = 1 : N
   ACFtab = ACFtabulate(ACF, k, c);
   ACCtab = ACFtabulate(ACC, k, c);

   if ~isempty(ACFfrac)
      ACFfractab = repmat(NaN, [size(ACFtab) size(ACFfrac, 4)]);
      for f = 1 : size(ACFfrac, 4)
         ACFfractab(:,:,f) = ACFtabulate(ACFfrac(:,:,:, f), k, c);
      end
   end

   if ~isempty(ACCfrac)
      ACCfractab = repmat(NaN, [size(ACCtab) size(ACCfrac, 4)]);
      for f = 1 : size(ACCfrac, 4)
         ACCfractab(:,:,f) = ACFtabulate(ACCfrac(:,:,:, f), k, c);
      end
   end

   if ~isempty(otherACF)
      for o = 1 : length(otherACF)
         ACFotherTab{o} = ACFtabulate(otherACF{o}, k, c);
      end
   end

   if ~isempty(otherACC)
      for o = 1 : length(otherACC)
         ACCotherTab{o} = ACFtabulate(otherACC{o}, k, c);
      end
   end

   for r = 1 : N

      if r < c
         CorrFlag = true;
      else
         CorrFlag = false;
      end

      subplot(N, N, (r - 1) * N + c)

      hold on

      if CorrFlag
         jack = plot(-k : k, fliplr(ACCtab(r, :)), 'b-', 'LineWidth', 3);
%          title(sprintf('cor(%s_t, %s_{t-k})', VariableNames{c}, VariableNames{r}))
         if ~isempty(VariableNames)
            title(sprintf('cor(%s_t, %s_{t-k})', VariableNames{r}, VariableNames{c}))
         end
         if ~isempty(ACCfrac)
            plotCI(fliplr(ACCtab(r, :)), flipud(squeeze(ACCfractab(r,:,:))),-k : k);
         end
         
         
         if ~isempty(otherACC)
            MultiColors = Colors4Plots;
            for o = 1 : length(otherACC)
               plot(-k : k,  fliplr(ACCotherTab{o}(r, :)), 'color', MultiColors{o+1});
            end
         end

         set(gca,'XTick', -k : tickSize : k)
         xlim([-k, k])
         ylim([-1, 1])
      else
         jack = plot(-k : k, ACFtab(r, :), 'b-', 'LineWidth', 3);
         if ~isempty(VariableNames)
            title(sprintf('cov(%s_t, %s_{t-k})', VariableNames{c}, VariableNames{r}))
         end
         if ~isempty(ACFfrac)
            plotCI(ACFtab(r, :), squeeze(ACFfractab(r,:,:)),-k : k);
         end

         if ~isempty(otherACF)
            MultiColors = Colors4Plots;
            for o = 1 : length(otherACF)
               plot(-k : k,  ACFotherTab{o}(r, :), 'color', MultiColors{o+1});
            end
         end

         if sameLimFlag
            ylim([miny, maxy])
         end
         set(gca,'XTick', -k : tickSize : k)
         xlim([-k, k])

         % make sure that zero is included in ylim
         ylm = ylim;
         if ylm(1) > 0
            ylim([0 ylm(2)]);
         end
         if ylm(2) < 0
            ylim([ylm(1) 0]);
         end
      end
      h(r,c) = gca;
   end
end

% draw origin
for r = 1 : N
   for c = 1 : N
      xlm = get(h(r,c), 'xlim');
      ylm = get(h(r,c), 'ylim');
      plot(h(r,c), xlm, [0 0], 'k:')
      plot(h(r,c), [0 0],  ylm, 'k:')
      set(h(r,c), 'xlim', xlm);
      set(h(r,c), 'ylim', ylm);
   end
end

if nargout == 0
   clear h
end