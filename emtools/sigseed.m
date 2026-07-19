function seed = sigseed(vals)
    % SIGSEED Deterministic hash seed from a vector of values.
    %   seed = SIGSEED(vals) returns a reproducible hash seed computed with a
    %   multiplicative rolling hash (Knuth's multiplicative constant, FNV
    %   offset basis as the initial accumulator), reduced mod a prime < 2^32
    %   each step. Input must be a nonempty vector of nonnegative integers.
    %
    %   Unlike keyHash, this produces identical seeds for the same input
    %   across MATLAB versions, so it is safe for archival reproducibility.
    %
    %   See also keyHash, RandStream, rng, rand, randi, randn

    assert(~isempty(vals), 'sigseed:emptyInput', ...
           'vals must be nonempty');
    assert(all(vals(:) >= 0 & vals(:) == floor(vals(:))), 'sigseed:badInput', ...
           'vals must be nonnegative integers');

    vals = transpose(uint64(vals(:)));  % row vector so the loop iterates element-by-element

    p = 4294967291;          % largest prime < 2^32
    m = uint64(2654435761);  % Knuth multiplicative constant
    h = uint64(2166136261);  % FNV offset basis (init only)
    for v = vals
        h = mod(h*m + v, p);   % mod every step keeps h < 2^32, so h*m never saturates uint64
    end
    seed = double(h);
end