function delta = spamaxabs(s)

delta = full(max(abs(s), [], 'all')); % enforcing full output to facilitate subsequent printing routines