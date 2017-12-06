if exist('wrap', 'var') && ~isempty(wrap)
   fprintf('\n');
   toc
   fprintf('\n');

   if isunix
       wrap = latexwrapper(wrap, 'compileDVI2PDF');
   else
       wrap = latexwrapper(wrap, 'close');
   end
   
   fprintf('Wrap file %s finished in %s\n', wrap.name, wrap.dir)
   fprintf('The current time is %s.\n', datestr(now))
end
