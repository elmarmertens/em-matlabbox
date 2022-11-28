function p = localresultsMCMC
% localresultsMCMC designates a folder where matlab outputs will be stored
% default: tmp folder in current working directory

p = fullfile(pwd, 'tmp');


if ~exist(p, 'dir')
    mkdir(p);
end
