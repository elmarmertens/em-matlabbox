function id = parid()
% PARID ...
%
%   ...

%% VERSION INFO
% AUTHOR    : Elmar Mertens
% $DATE     : 21-Dec-2017 09:25:27 $
% $Revision : 1.00 $
% DEVELOPED : 9.2.0.556344 (R2017a)
% FILENAME  : parid.m

try 
    task = getCurrentTask;
    % Please see the "Parallel Computing Toolbox" section of the matlab preferences to see whether a pool will be automatically created. 
    % To open these preferences run "preferences('Parallel Computing Toolbox')" from the Matlab command line.

catch poolME
    % one typical cause of problems could be limited availability of licenses
    warning(poolME.identifier, 'There was a problem obtaining a pool of parallel workers; trying to work without one.\n The error message was:\t %s', ...
        poolME.message)
    task = []; % gcp('nocreate'); % If no pool can be created, do not create new one.
end

if isempty(task)
    id = 1;
else
    id = get(task, 'ID');
end
