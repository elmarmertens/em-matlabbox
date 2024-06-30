function is_m1_mac = isarm64()
% ISARM64 ... 
%  
%   ... 

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 16-Jan-2024 16:35:34 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 23.2.0.2459199 (R2023b) Update 5 
% FILENAME  : isarm64.m 

if ismac()
    [~,result] = system('uname -v');
    is_m1_mac = any(strfind(result,'ARM64'));
else
    is_m1_mac = false;
end