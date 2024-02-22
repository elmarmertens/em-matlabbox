function [KSC, KSCt, logy2offset] = getKSC10values(T, Nsv)
% GETKSC10VALUES returns coefficients from Omori, Chib, Shephard, Nakajima 10-point extension
% to the original 7-point mixture of Kim, Shephard and Chib to approximate a chi2 variable
%
% USAGE: [KSC, KSCt, logy2offset] = getKSC10values(T, Nsv) returns in addition, the structure KSCt, with the same fieldnames as KSC
%
% where KSC returns a structure with elements mean, vol, var, pdf, and cdf
% KSCt returns are corresponding structure but with fields "blown up" to dimension T x Nsv x 10
% and logy2offset is the offset c used to compute log(y^2 + c)
%

%   Coded by  Elmar Mertens, em@elmarmertens.com


if nargin < 2
    Nsv = 1;
end

KSC.mean    = [1.92677 1.34744 0.73504 0.02266 -0.85173 -1.97278 -3.46788 -5.55246 -8.68384 -14.65000];
KSC.var     = [0.11265 0.17788 0.26768 0.40611 0.62699 0.98583 1.57469 2.54498 4.16591 7.33342];
KSC.vol     = sqrt(KSC.var);
KSC.pdf     = [0.00609 0.04775 0.13057 0.20674 0.22715 0.18842 0.12047 0.05591 0.01575 0.00115];
KSC.cdf     = cumsum(KSC.pdf);

% blowup to cover time dimension
if nargout > 1
    fn = fieldnames(KSC);
    for f = 1 : length(fn)
        KSCt.(fn{f}) = repmat(permute(KSC.(fn{f}), [1 3 2]), [Nsv, T, 1]);
    end
end

logy2offset = 0.0001; % see OCSN page 436, this is the offset c used to compute log(y^2 + c)
