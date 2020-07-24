function [PAIdraws, cholSIGMAdraws] = BVARdraws(PAIhat, sqrtXX, sqrtSSR, dof, Ndraws, rndStream)
% BVARDRAWS ... 
%  
%   ... 

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 03-Jun-2020 22:08:32 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 9.8.0.1380330 (R2020a) Update 2 
% FILENAME  : BVARdraws.m 


if nargin < 5
    Ndraws = 1;
end
if nargin < 6
    rndStream = getDefaultStream;
end

[K, N] = size(PAIhat);

cholSIGMAdraws = iwishcholdraw(sqrtSSR, dof, Ndraws, rndStream); 

PAIdraws   = NaN(K,N,Ndraws);
vecPAI     = PAIhat(:);
invsqrtXX  = eye(K) / sqrtXX;

zdraws     = randn(rndStream, N * K, Ndraws);

for n = 1 : Ndraws
    sqrtOMEGApai    = kron(cholSIGMAdraws(:,:,n), invsqrtXX'); % note: not lower triangular due to transpose on invsqrtXX
    thisPAI         = vecPAI + sqrtOMEGApai * zdraws(:,n);
    PAIdraws(:,:,n) = reshape(thisPAI, K, N);
end
