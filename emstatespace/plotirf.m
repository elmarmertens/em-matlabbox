function h = plotirf(irf, ynames, wnames, irf2, ticksize, plotname, linetype)
% function h = plotirf(irf, ynames, wnames, irf2, ticksize, plotname)
% plots irf and labels rows with ynames and columns with wnames
%
% if irf2 is specified, it will be superimposed on the plots next to irf
% (irf2 needs to be of same dimension as irf)
%
% ticksize: specify steps of xtick-labels, default: 4, use auto-setting if empty
%
% h returns ny x nw array of handels for each irf's subplot

%   Coded by  Elmar Mertens, em@elmarmertens.com

% Elmar Mertens
% www.elmarmertens.ch

narginchk(1, 7)

[ny, nw, lags] = size(irf);
lags = lags - 1; % time starts at zeros

if nargin < 2
   ynames = [];
end

if nargin < 3
   wnames = [];
end

if nargin < 4
   irf2 = [];
elseif ~isempty(irf2) 
   tmp1 = size(irf2);
   tmp2 = size(irf);
   if ~isequal(tmp1(1:3),tmp2(1:3))
       error('irf2 and irf need to share the same dimension (at least one to three)')
   end
end

if nargin < 5 || isempty(ticksize)
   ticksize = 4; % typically IRF's are quarterly
   % ticksize = max(2,floor(lags / 5 / 4) * 4); % alternative setting
end

if nargin < 6
    %    if ~isempty(inputname(1))
    %       plotname = inputname(1);
    %    else
    %       plotname = 'IRF';
    %    end
    plotname = [];
end

if nargin < 7
    linetype = {'b-', 'LineWidth', 2};
end

timeAxis = 0 : lags;
if nargout > 0
   h = zeros(ny, nw);
end

% clf

% loop over rows and columns
for r = 1 : ny
   for c = 1 : nw

      subplot(ny, nw, (r - 1) * nw + c);
      hold on

      plot(timeAxis', squeeze(irf(r, c, :)), linetype{:});
      if ~isempty(irf2)
          colors = Colors4Plots;
          colors = colors(2:end);
          for n = 1 : size(irf2,4)
              plot(timeAxis', squeeze(irf2(r, c, :,n)), '-', 'color', colors{n}, 'LineWidth', 2);
              grid on
          end
      end
      if nargout > 0
         h(r, c) = gca;
      end

      if c == 1 && ~isempty(ynames)
         ylabel(ynames(r));
      end
      if r == 1 && ~isempty(wnames)
         title(wnames(c));
      end

      % draw zero line
      plot(timeAxis([1 end]), [0 0], 'k:')
      xlim([0 lags])
      if ~isempty(ticksize)
         set(gca,'XTick', timeAxis(1:ticksize:end))
      end
   end
   
end

if ~isempty(plotname)
    set(gcf, 'name', plotname)
end
