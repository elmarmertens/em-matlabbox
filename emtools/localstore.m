function p = localstore
% LOCALSTORE designates a folder where all matlab outputs will be stored
% default: tmp folder in current working directory

p = fullfile(pwd, 'tmp');


if ~exist(p, 'dir')
    mkdir(p);
end
