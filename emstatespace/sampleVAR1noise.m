function [h, hbar, hshock, htilde] = sampleVAR1noise(obs, Nsv, T, rho, hVCVsqrt, Eh0, sqrtVh0, noisevol, rndStream)

zerosNsv  = zeros(Nsv);
Insv      = eye(Nsv);
A         = [diag(rho) zerosNsv; zerosNsv Insv];
B         = [hVCVsqrt; zerosNsv];
C         = [Insv Insv];
sqrtR = zeros(Nsv,Nsv,T);
for n = 1 : Nsv
    sqrtR(n,n,:) = noisevol(n,:);
end


x0           = [zeros(Nsv, 1); Eh0];
sqrtVhtilde  = eye(Nsv); % Note: need fixed prior, not dependent on estimated rhos (alt: use prior rho)
sqrtVx0      = [sqrtVhtilde, zerosNsv; zerosNsv sqrtVh0];
[H, Hshock, H0] = a2b2c2DisturbanceSmoothingSampler1draw(A, B, C, obs, x0, sqrtVx0, ...
    sqrtR, rndStream); 

h      = H(1:Nsv,:) + H(Nsv+1:end,:); % C * H
hbar   = H0(Nsv+1:end);
htilde = cat(2, H0(1:Nsv,:), H(1:Nsv,:)); % demeaned component, including lagged value
hshock = Hshock(1:Nsv,:);
