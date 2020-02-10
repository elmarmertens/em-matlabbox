function wrapopen(wrap)
% opens compiled wrap file (if it exists)

if isempty(wrap) || ~(wrap.id == 0)
    return
end

if ismac
    system(sprintf('open %s.pdf', fullfile(wrap.dir, wrap.name)))
end