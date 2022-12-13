if exist('wrap', 'var') && ~isempty(wrap)
    fprintf('\n');
    toc
    fprintf('\n');

    %     if ispc
    %         wrap = latexwrapper(wrap, 'compilepdf');
    %     else
    %         wrap = latexwrapper(wrap, 'compileDVI2PDF');
    %     end

    wrap = latexwrapper(wrap, 'compileDVI2PDF');

    % if compilation fails on your system, simply replace the line above with the following:
    % wrap = latexwrapper(wrap, 'close');
    
    fprintf('Wrap file %s finished in %s\n', wrap.name, wrap.dir)
    fprintf('The current time is %s.\n', datestr(now))
end
