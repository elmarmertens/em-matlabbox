function lambda = ploteigenvalues(A, B)
% lambda = function ploteigenvalues(A,B)
% function ploteigenvalues(lambda)
% 
% plot lambdas over the unit circle, Inf eigenvalues are omitted

% Elmar Mertens
% www.elmarmertens.ch

error(nargchk(1,2,nargin))
if nargin == 2
   lambda = eig(B, A);
else 
   lambda = A;
end

figure
% set(gca, 'fontsize', 16)
hold on

%first unit circle
theta    = 0:.1:360;
[X, Y]   = pol2cart(theta, 1);
plot(X,Y, 'r-')

%now the lambdas
plot(real(lambda), imag(lambda), 'bx')
plot(real(lambda), imag(lambda), 'bo')

% the next lines would mark instable roots with a diamond:
% outside = abs(lambda) > 1;
% if any(outside)
%     plot(real(lambda(outside)), imag(lambda(outside)), 'ks', ...
%         'MarkerEdgeColor','b',...
%         'MarkerFaceColor','b',...
%         'MarkerSize',10)
% end
    
hold off
maxx = max(abs([xlim, ylim]));

xlim([-maxx maxx])
ylim([-maxx maxx])

ylabel('imag(\lambda)')
xlabel('real(\lambda)')
plotOrigin
axis square
if nargin == 1
   ttl = '\lambda';
   set(gcf, 'name', 'eigenvalues')
else
   ttl = '|{\bf A} \lambda - {\bf B}| = 0';
   set(gcf, 'name', 'eig(B, A)')
end
if any(abs(lambda) == Inf)
   ttl = sprintf('%s\n(%d infinite eigenvalues)', ttl, sum(abs(lambda) == Inf));
end
title(ttl)