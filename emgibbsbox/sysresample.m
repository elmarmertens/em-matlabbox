function ndx = sysresample(pdf, N, thisudraw, iscdf)
% RANDNDXMULTINOMIAL ... 
%  
%   ndx = systesample(pdf, N, iscdf)
%   ... 

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 18-Dec-2013 14:50:30 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 7.14.0.739 (R2012a) 
% FILENAME  : randndxmultinomial.m 

if nargin < 4
    iscdf = false;
end


if iscdf
    cdf = pdf(:);
else
    cdf = cumsum(pdf(:));
end


udraws = (thisudraw + (1:N) - 1) / N;

%% exploit that udraw and cdf are sorted
ndx   = NaN(1, N);
this = 1;
for n = 1 : N    
    
    while cdf(this) < udraws(n)
        this = this + 1;
    end
        
    ndx(n) = this; 
end

%% check
% ndx0    = sum(bsxfun(@gt, udraws, cdf),1) + 1;
% if ~isequal(ndx, ndx0)
%     warning('something off')
% end
% bsxfun seems faster
% ndx2   = sum(udraw > cdf, 1) + 1;
% checkdiff(ndx, ndx2);
