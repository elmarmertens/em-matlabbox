if true || isunix % && ~isdesktop
    clear wrap
    
    wrap.dir    = fullfile(localtemp, 'foo');
    if ~exist('titlename', 'var')
        titlename   = mfilenamecaller;
    end
    if isempty(titlename)
        wrap.title  = sprintf('\\titlecaveat{Dummy wrapper}');
        wrap.name   = 'foo';
    else
        latextitle = strrep(titlename, '_', ' ');
        % insert a linebreak if title too long
        if length(latextitle) > 50
            dashes = strfind(latextitle, '-');
            ndx = find(dashes <= 50, 1, 'last');
            latextitle = [latextitle(1:dashes(ndx)), '\ldots\linebreak ', ...
                latextitle(dashes(ndx)+1:end)];
        end
        wrap.title = sprintf('\\titlecaveat{%s}', latextitle);
        
        % remove any LaTeX line breaks included in titlename
        % titlename = strrep(titlename, '\\', '');
        
        [jim, titlename] = fileparts(titlename);
        wrap.name    = titlename;
    end
    if exist('doCharts', 'var')
        wrap.doDcolColors = true;
    end
    % wrap.name = strcat(wrap.name, datestr(now, 30));
    wrap = latexwrapper(wrap, 'start');
    % wrap = diary2wrap(wrap, [], true);
    tic
else
    wrap = [];
end
