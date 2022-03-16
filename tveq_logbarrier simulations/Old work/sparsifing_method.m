
clear
close all force

% path(path, './Optimization');
% path(path, './reconstruction');
% path(path, './Measurements');
% path(path, './Data');
% path(path, './Usefull stuff');

addpath('C:\Users\Dvir\Desktop\Doppler Frequencies Mapping')
addpath('C:\Users\Dvir\Desktop\Doppler Frequencies Mapping\Compressive_Sensing_Tools\l1magic\Optimization')
addpath 'C:\Users\Dvir\Desktop\Doppler Frequencies Mapping\Usefull stuff'
addpath 'C:\Users\Dvir\Desktop\DMD Settings and Control\Compressive_Sensing_Tools\NESTA_v1.1\Analysis'


addpath utils/
addpath utils/qr/
addpath utils/utils_Wavelet

% No. of reweighting iterations (Optional).
rwt = 5;       
pow = 5;
sz = 2^pow;

N = sz^2; 

% No. of measurements
M = round(N/8); 

xs = linspace(-1, 1, sz);
ys = xs;

[xx, yy] = meshgrid(xs, ys);

w = 0.4;
sig_rad = 0.4;

t = 1;

% Refferenced signal.
signal_type = 1;

switch signal_type
    case 1
        % Gaussian signal with incertainty of sqrt(w).
        r = exp(-((xx.^2) + (yy.^2)) ./ w).*exp(-1i.*pi*t);
        
    case 2
        % Circular signal radius sig_rad.
        r = circle_draw(sz, sig_rad);
end

Ir = abs(r).^2; 

% Refferenced signal (has no usge now).
s = r;
Is = abs(s).^2; 

photons = 2*sum(Ir, 'all');

% Laser wavelenght.
lambda = 532e-9;

% Piazo engeen doppler frequencies membren.
membrane_type = 1;
membrane_scale = 1e-8;

switch membrane_type
    case 0 
        % Gaussian membrane.
        vD = exp(-(xx.^2 + yy.^2)./0.4) * membrane_scale;
        vD(xx.^2 + yy.^2 > 0.65^2)= 0;

    case 1
        % Half circle membrane.
        R = 0.7;
        vD = half_circle_draw(sz, R) * membrane_scale;
        
    case 2
        % Membrane shape.
        vD = membren(sz).* 10^-8;
        
    case 3
        % Advanced membrane (similar to real one).
        vD = advanced_membrane(sz) * membrane_scale;

    case 4
        % Advanced membrane (Different version).
        vD = advanced_membrane2(sz) * membrane_scale;

    case 5 
        % Advanced membrane with control of:
        angs = 5; % N. of angels.
        rads = 3; % N. of radiuses.
        
        vD = advanced_membrane_set(sz, angs, rads) * membrane_scale;
        
    case 6
        % Parabul membrane.
         vD = (xx.^2 + yy.^2) * membrane_scale;
         vD(xx.^2 + yy.^2 > 0.65^2)= 0;

   otherwise
        disp('Invalid option for "membrane" variable');
end
pcolor(xs, ys, vD);

% Doppler frequencies.
fD = 2.*vD./lambda; 

% PHYSICAL PARAMETERS (Raman-Nath Scattering (lamda*L / d^2 << 1): 

% Beams wave lenghts, in [m].
lambda = 532e-9; 

% Magnitude of the laser wavevector, in [1/m].
k0 = 2*pi / lambda; 

% In case of same intens. The intensity of each beam in [mW/cm^2].
I = 0.3; 

% Refractive index.
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

% After liquid crystal 
I1 = (abs(besselj(0, rho) + 1i.*besselj(1, rho).*exp(-1i.*psi)).^2).*Ir;
I0 = (abs(besselj(-1, rho) + 1i.*besselj(0, rho).*exp(-1i.*psi)).^2).*Ir;

photons_after_liquid = sum(I1, 'all') + sum(I0, 'all');

balance = I1 - I0;

flat_balance = reshape(balance.',1,[])';
haar = generate_haar(sz^2);
sprs = haar * flat_balance;

sps_energy = sum(sprs, 'all');
balance_energy = sum(balance, 'all');

[val, ind] = sort(abs(sprs),'descend');
    ind_pos = ind(val>1e-4);
gamma_orig = ind_pos(1:min(length(ind_pos),M-1));

% x_gues = zeros(N,1);
%  xAbs = sort(abs(x_gues),'descend');
%     cutoff = xAbs( round( .02*N) );
%     u = 1./( abs(x_gues) + cutoff );
%     U = spdiags( u, 0, N, N );
%     Ut = U;
%     normU = max(u);


