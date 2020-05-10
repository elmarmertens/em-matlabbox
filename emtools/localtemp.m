function p = localtemp
% LOCALSTORE designates a folder where all matlab store temp files
% default: tmp folder in current working directory

p = pwd;


if ~exist(p, 'dir')
    mkdir(p);
end