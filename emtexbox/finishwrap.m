if exist('wrap', 'var') && ~isempty(wrap)
    fprintf('\n');
    toc
    fprintf('\n');

        if ispc
            wrap = latexwrapper(wrap, 'close');
            latexwrapper(wrap, 'pdflatex');
            latexwrapper(wrap, 'pdflatex');
        elseif ismac
            wrap = latexwrapper(wrap, 'compileDVI2PDF');
        else
            wrap = latexwrapper(wrap, 'close');
        end

    % if compilation fails on your system, simply replace the line above with the following:
    % wrap = latexwrapper(wrap, 'close');
    
    fprintf('Wrap file %s finished in %s\n', wrap.name, wrap.dir)
    fprintf('The current time is %s.\n', datetime)
end
