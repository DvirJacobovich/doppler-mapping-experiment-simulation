% A_noiselet.m
%
% Takes noiselet measurements on the specified set.
%
% Written by: Justin Romberg, Georgia Tech, jrom@ece.gatech.edu
% Created: May 2007
%

function y = A_noiselet_3D(x, OMEGA, m, n, l)

N = length(x);

x_3d = reshape(x, m * n, l);
w = realnoiselet(x_3d) / sqrt(N);

y = w(OMEGA,:);
y = y(:);


