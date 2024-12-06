function [h, hbar, hshock, htilde] = sampleVAR1noisePBS(obs, Nsv, T, rho, hVCVsqrt, Eh0, sqrtVh0, noisevol, rndStream)


NsvT    = Nsv * T;
Nstates = Nsv * (T + 2); % plus initial values and shocks

% construct AA
rows1   = 1:Nstates;
cols1   = 1:Nstates;
values1 = ones(Nstates,1);

rowsRHO   = Nsv + (1 : NsvT);
colsRHO   = 1 : NsvT;
valuesRHO = -1 .* rho .* ones(1,T); % array expansion is faster than repmat(rho, 1, T)

rowsCONST   = rowsRHO;
colsCONST   = NsvT + Nsv + repmat((1 : Nsv)', 1, T); 
valuesCONST = -1 - valuesRHO; % this yields "-(1-rho) = -1 + rho")

rows   = cat(1, rows1(:), rowsRHO(:), rowsCONST(:));
cols   = cat(1, cols1(:), colsRHO(:), colsCONST(:));
values = cat(1, values1(:), valuesRHO(:), valuesCONST(:));
AA     = sparse(rows, cols, values, Nstates, Nstates);

% construct invBB
Insv        = eye(Nsv);
rows0       = 1:Nsv;
cols0       = rows0;
values0     = ones(Nsv,1);
rowsCONST   = Nsv + NsvT + repmat((1:Nsv)', 1, Nsv);
colsCONST   = Nsv + NsvT + repmat(1:Nsv, Nsv, 1);
valuesCONST = Insv / sqrtVh0;
rowsTilde   = Nsv + repmat((1:Nsv)', 1, Nsv) + Nsv .* permute(0:T-1, [1 3 2]);
colsTilde   = Nsv + repmat(1:Nsv, Nsv, 1) + Nsv .* permute(0:T-1, [3 1 2]);
valuesTilde = repmat(Insv / hVCVsqrt, 1, 1, T);
rows        = cat(1, rows0(:), rowsTilde(:), rowsCONST(:));
cols        = cat(1, cols0(:), colsTilde(:), colsCONST(:));
values      = cat(1, values0(:), valuesTilde(:), valuesCONST(:));
invBB       = sparse(rows, cols, values, Nstates, Nstates);

% construct CC
rows   = 1:NsvT;
cols   = Nsv + (1:NsvT);
values = ones(NsvT,1);
CC     = sparse(rows, cols, values, NsvT, Nstates);

% construct DD
invDD = spdiags(1 ./ noisevol(:),0,NsvT,NsvT);

% construct X0
XX0   = sparse(NsvT+Nsv+(1:Nsv), 1, Eh0, Nstates, 1);

%% means and innovations
AAtilde            = invBB * AA;
XX0tilde           = invBB * XX0;

CCtilde            = invDD * CC;
Ytilde             = invDD * obs(:);

P                  = AAtilde' * AAtilde + (CCtilde' * CCtilde);
sqrtP              = chol(P, 'lower');


%   [sqrtP, flag]       = chol(P, 'lower');
% if flag > 0
%     error('P not posdf, using QR instead')
% end

sqrtPXhat    = sqrtP \ (AAtilde' * XX0tilde + CCtilde' * Ytilde);
Zdraw        = randn(rndStream, Nstates, 1);
Xdraw        = transpose(sqrtP) \ (sqrtPXhat + Zdraw);
XshockDraw   = AA * Xdraw - XX0;

%% map into outputs
hndx     = Nsv + (1:NsvT);
h        = reshape(Xdraw(hndx), Nsv, T);

hbar     = Xdraw(NsvT + Nsv + (1:Nsv));

% htilde (including lag)
htilde   = reshape(Xdraw(1:Nsv+NsvT), Nsv, T+1) - hbar;

% hshock
hshock   = reshape(XshockDraw(hndx), Nsv, T);

end