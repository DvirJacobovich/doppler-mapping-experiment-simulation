% At_noiselet.m
%
% Adjoint of A_noiselet.m
%
% Written by: Justin Romberg, Georgia Tech, jrom@ece.gatech.edu
% Created: May 2007
%

function v = At_noiselet_3D(y, OMEGA, m, n, l, M)

y_mat = reshape(y,M,l);

y_full = zeros(m*n,l);
y_full(OMEGA,:) = y_mat;


v = realnoiselet(y_full/sqrt(M));
v = v(:);


