function [Xdraw, QQ, RR1, Xhat] = abcYXprecisionsampler(Y,XX0,AA,invBB,CC,rndStream,QQ,RR1)
% YXprecisionsampler bare-bones precision sampler for vectorized inputs.
%

%% VERSION INFO
% AUTHOR    : Elmar Mertens

% get dimensions
if nargin < 7
    QQ  = [];
    RR1 = [];
end


%% CC and prepare Arows and Brows

if isempty(QQ)
    % perform QR
    [QQ,RR]   = qr(CC');
    [N1, N2]  = size(CC);
    N2        = N2 - N1;
    RR1       = RR(1:N1,1:N1)';
else
    N1        = size(RR1,1);
    N2        = size(QQ,1) - N1;
end

QQ1       = QQ(:,1:N1)';
QQ2       = QQ(:,N1+1:end)';


%% means and innovations
EX        = AA \ XX0;
EY        = CC * EX;

X1tilde   = RR1 \ (Y - EY);

QQX1tilde = QQ1' * X1tilde;

%% precision-based sampler
AAtilde       = invBB * AA;
AAtildeQQX1   = AAtilde * QQX1tilde;
AAtildeQQ2    = AAtilde * QQ2';
% invQSIG21     = AAtildeQQ2' * AAtildeQQX1;
invQSIG22     = transpose(AAtildeQQ2) * AAtildeQQ2;
cholinvQSIG22 = chol(invQSIG22, 'lower');

X2hat         = - cholinvQSIG22 \ (AAtildeQQ2' * AAtildeQQX1);

Z2draw        = randn(rndStream, N2, 1) + X2hat;
X2draw        = cholinvQSIG22' \ Z2draw;
Xdraw         = EX + QQX1tilde + QQ2' * X2draw;

if nargout > 3
    Xhat       = EX + QQ1' * X1tilde + QQ2' * (cholinvQSIG22' \ X2hat);
end
