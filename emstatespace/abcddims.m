function abcd = abcddims(abcd)
% ABCDDIMS
% sets dimensions of abcd
% USAGE: abcd = abcddims(abcd)

%   Coded by  Elmar Mertens, em@elmarmertens.com


if ~isfield(abcd, 'C') && ~isfield(abcd, 'D') && isfield(abcd, 'Cx')
   abcd.C = abcd.Cx * abcd.A;
   abcd.D = abcd.Cx * abcd.B;
end

error(abcdchk(abcd.A,abcd.B,abcd.C,abcd.D));
% [abcd.chkmsg,abcd.A,abcd.B,abcd.C,abcd.D] = abcdchk(abcd.A,abcd.B,abcd.C,abcd.D);

abcd.Ny           = size(abcd.C, 1);
[abcd.Nx abcd.Nw] = size(abcd.B);

if isfield(abcd, 'Xnames') && length(abcd.Xnames) ~= abcd.Nx
   warning('em:msg', '# of Xnames (%d) not consistent with Nx=%d', length(abcd.Xnames), abcd.Nx);
end
if isfield(abcd, 'Wnames') && length(abcd.Wnames) ~= abcd.Nw
   warning('em:msg', '# of Wnames (%d) not consistent with Nw=%d', length(abcd.Wnames), abcd.Nw);
end
if isfield(abcd, 'Ynames') && length(abcd.Ynames) ~= abcd.Ny
   warning('em:msg', '# of Ynames (%d) not consistent with Ny=%d', length(abcd.Ynames), abcd.Ny);
end

%% prune old results
if isfield(abcd, 'EY')
   abcd = rmfield(abcd, 'EY');
end
if isfield(abcd, 'EX')
   abcd = rmfield(abcd, 'EX');
end
if isfield(abcd, 'IRFy')
   abcd = rmfield(abcd, 'IRFy');
end
if isfield(abcd, 'IRFx')
   abcd = rmfield(abcd, 'IRFx');
end
if isfield(abcd, 'GammaX')
   abcd = rmfield(abcd, 'GammaX');
end
if isfield(abcd, 'GammaY')
   abcd = rmfield(abcd, 'GammaY');
end
if isfield(abcd, 'GammaXcond')
   abcd = rmfield(abcd, 'GammaXcond');
end
if isfield(abcd, 'GammaYcond')
   abcd = rmfield(abcd, 'GammaYcond');
end
