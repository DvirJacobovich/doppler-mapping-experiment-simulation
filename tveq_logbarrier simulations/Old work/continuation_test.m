
clear
close all force

path(path, './Optimization');
path(path, './Measurements');
path(path, './Data');
addpath('C:\Users\Dvir\Desktop\Doppler Frequencies Mapping')

addpath utils/
addpath utils/qr/
addpath utils/utils_Wavelet

% reweighted setup
rwt = 5;        % number of reweighting iterations
pow = 6;
sz = 2^pow;
N = sz^2; 

% No. of measurements
M = round(N/8); 

% Sparsity level
T = round(M/3);    

xs = linspace(-1, 1, sz);
ys = xs;

[xx, yy] = meshgrid(xs, ys);
w = 0.4;
t = 1;

% Refferenced signal.

signal = 1;

switch signal
    case 1
        % Gaussian signal multiplied by spatial-phase.
        r = exp(-((xx.^2) + (yy.^2)) ./ w).*exp(-1i.*pi*t);
        
    case 2
        % Circular signal 0.4 radius.
        r = circle_draw(sz, 0);
end

Ir = abs(r).^2; 

% Refferenced signal (has no usge now).
s = r;
Is = abs(s).^2; 

photons = 2*sum(Ir, 'all');

% Piazo engeen doppler frequencies membren.
lambda = 532e-9;

membrane = 3;

switch membrane
    case 0 
        % Gaussian membrane.
        vD = exp(-(xx.^2 + yy.^2)./0.4).*10^-8;
        vD(xx.^2 + yy.^2 > 0.65^2)= 0;

    case 1
        % Circular membrane. for half circule - 1, for full circule - 0.
        vD = circule_draw2(sz, 1) * 10^-8;
        
    case 2
        % Membrane shape.
        vD = membren(sz);
        
    case 3
        % Advanced membrane (similar to real one).
        vD = advanced_membrane(sz) .* 10^-8;

    case 4
        % Advanced membrane (Different version).
        vD = advanced_membrane2(sz) .* 10^-8;

    case 5
        % Parabul membrane.
         vD = (xx.^2 + yy.^2).*10^-8;
         vD(xx.^2 + yy.^2 > 0.65^2)= 0;

   otherwise
        disp('Invalid option for "membrane" variable');
end
pcolor(xs, ys, vD);

% Doppler frequencies.
fD = 2.*vD./lambda; 


% Physical parameters, Raman-Nath Scattering (lamda*L / d^2 << 1).

% Beams wave lenghts, in [m].
lambda = 532e-9; 

% Magnitude of the laser wavevector, in [1/m].
k0 = 2*pi / lambda; 

% In case of same intens. The intensity of each beam in mW/cm^2.
I = 0.3; 

n2 = -1.8*10^-4; 

% Fringe spacing.
d = 50*10^-6; 

% Interaction region thickness.
L = 9*10^-6; 

ld = 12*10^-6; 

% Grating index wavevector.
kg = 2*pi/(110*10^-6); 

% Liquid crystal relaxation time.
tau = 72*10^-3; 

rho = (2*k0*n2.*d ./ (((1 + ld^2*kg^2)^2 + (fD.*tau)^2).^0.5)).*(Is.*Ir).^0.5;
psi = atan(fD.*tau ./ (1 + ld^2*kg^2)^2);

% after liquid crystal 
I1 = (abs(besselj(0, rho) + 1i.*besselj(1, rho).*exp(-1i.*psi)).^2).*Ir;
I0 = (abs(besselj(-1, rho) + 1i.*besselj(0, rho).*exp(-1i.*psi)).^2).*Ir;

photons_after_liquid = sum(I1, 'all') + sum(I0, 'all');

balance = I1 - I0;
figure, image(balance,'CDatamapping','scaled'),set(gca,'FontSize',20);
title('Original balanced signal');

% We start with one white balance measurement.
white_bal_measre = sum(balance, 'all');

permutation1 = randperm(N, M);
permutation0 = randperm(N, M);

% Two different masks matrices, each one is positive real dragon wavelet.
masks1 = double(masksMat(N, ...
    M, permutation1)');

masks0 = double(masksMat(N, ...
    M, permutation0)');

b1 = masks1 * reshape(I1.',1,[])';
b0 = masks0 * reshape(I0.',1,[])';

% First balanced projection.
b = b1 - b0; 

y1 = masks0 * reshape(I1.',1,[])';
y0 = masks1 * reshape(I0.',1,[])';

% Second balanced projection.
y = y1 - y0; 

% Subtracting the White Balanced measurment.
measured = y + b - white_bal_measre; 
A = masks1 + masks0 + 1;

At = A';

% Final measurement matrix pointer function.
A_f = @(x)A * x;
At_f = @(x)At * x;

% Initial gues.
x0 = At * measured; 

tvI = sum(sum(sqrt([diff(balance,1,2) zeros(sz,1)].^2 + [diff(balance,1,1); zeros(1,sz)].^2 )));
fprintf('Original TV = %.3f\n', tvI);

time0 = clock;

% Algorithm terminates when the duality gap <= duality_gap_barrier.
duality_gap_barrier = 1e-4; 

% Tolerance for SYMMLQ.
slqtol = 1e-9;

for i = 1:3
    xp =  tveq_logbarrier(x0, A_f, At_f, measured, duality_gap_barrier, 10, slqtol*10^(-i), 200); 
    Ip = reshape(xp, sz, sz)';
    figure, image(Ip,'CDatamapping','scaled'),set(gca,'FontSize',20);
    title('rec balanced Pos/Neg fix. Per 12.5%', i);
    x0 = xp;
end

fprintf('Total elapsed time = %f secs\n\n', etime(clock,time0));

% Error calculations.

% 1. General relative (per pixel) error.
D = abs(Ip - balance).^2;
relative_MSE_err = sum(D(:)) / numel(Ip);


% 2. Only relevant area error.
if signal == 1 % Gaussian signal with standard deviation w.
    center = xx.^2 + yy.^2 < w;
else 
    center = xx.^2 + yy.^2 < 0.4; % Circular signal with 0.4 radius.
end

rel_Ip = Ip .* center;
rel_balance = balance .* center;

% Number of elements to be taken.
flat_cen = reshape(center.',1,[])';
num = size(flat_cen(flat_cen ~= 0));
num = num(1);

D_rel = abs(rel_Ip - rel_balance).^2;
relative_rel = sum(D_rel(:)) / num;

xi = 8*pi*tau*besselj(0, 2*k0*d*n2.*Ir)*besselj(1, 2*k0*d*n2.*Ir);

del = Ip ./ Ir;
dni = del ./ xi;

% figure, image(dni,'CDatamapping','scaled'),set(gca,'FontSize',20);
% title('dni');

                                                                                                                         
