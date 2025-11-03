function hl = plotPITcdfstairs(thesepits, thesehalfbands, theselabels, doLegendBands)

if nargin < 4 || isempty(doLegendBands)
    doLegendBands = false;
end

% fontsize   = 18;
Nmodels    = length(theselabels);
% rvec    = 0 : 0.001 : 1;
legHandles = [];
legLabels  = cell(0);
% thisfig    = figure;
% set(gca, 'FontSize', fontsize);
hold on;
for m = 1 : Nmodels
    thispit   = thesepits(:,m);
    if ~all(isnan(thispit))
        [f,x]     = ecdf(thispit);
        this      = stairs(x,f, 'LineWidth', 2, 'Color', colors4plots(m));
        thishalfband = thesehalfbands(m);
        legHandles = cat(1, legHandles, this);
        legLabels  = cat(1, legLabels, theselabels(m));
        if contains(theselabels(m), 'VAR')
            this.LineStyle = '-.';
        end
        % get bands
        [~,~,~,~,KS95cv,Tpit,rvec] = PITtest(thispit); %#ok<ASGLU>
        bplus  = plot(rvec,rvec + thishalfband,'Color',this.Color,'LineStyle','-','LineWidth',1);
        bminus = plot(rvec,rvec - thishalfband,'Color',this.Color,'LineStyle','-','LineWidth',1);
        if contains(theselabels(m), 'VAR')
            bplus.LineStyle = '-.';
            bminus.LineStyle = '-.';
        end
        if doLegendBands
            legHandles = cat(1, legHandles, bplus);
            legLabels  = cat(1, legLabels, sprintf('Confidence band (%s)', theselabels{m}));
        end
    end
end
plotfortyfive('k:', 'linewidth', 2);
xlim([0 1]); ylim([0 1]);
% legHandles = legHandles([1 4 2]);
% legLabels  = legLabels([1 3 2]);

hl = legend(legHandles, legLabels, 'location', 'northwest');
% if ~isempty(thistitle)
%     title(thistitle)
% end

end % function plotPITcdfstairs