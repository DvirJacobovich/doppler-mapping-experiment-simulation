% Dvir Jacobovich 2021 - Pr John Howell lab Hebrew University Of Jerusalem

clear
close all force

addpath('..\')

addpath utils/
addpath utils/qr/
addpath utils/utils_Wavelet

% No. of reweighting iterations (Optional).
pow = 7;
sz = 2^pow;

N = sz^2; 

% No. of measurements
M = round(N); 

xs = linspace(-1, 1, sz);
ys = xs;

[xx, yy] = meshgrid(xs, ys);

Amp = 50*1e-2; % 50mV
w = 0.6;
sig_rad = 0.6;

t = 1;

% Refferenced signal.
signal_type = 1;

switch signal_type
    case 1
        % Gaussian signal with incertainty of sqrt(w).
        r = Amp*exp(-((xx.^2) + (yy.^2)) ./ w).*exp(-1i.*pi*t);
        
    case 2 % Uniformed intensity
        % Circular signal radius sig_rad.
        r = Amp*circle_draw(sz, sig_rad);
end

Ir = abs(r).^2; 

% Refferenced signal (has no usge now).
s = r;
Is = abs(s).^2; 

photons = 2*sum(Ir, 'all');

% Beams wave lenghts, in [m].
lambda = 532e-9;

% Piazo engeen doppler frequencies membren.
piezo_type = 7;
velocity_scale = 1e-8;

switch piezo_type
    case 0 
        % Gaussian membrane.
        vD = exp(-(xx.^2 + yy.^2)./0.4) * velocity_scale;
        vD(xx.^2 + yy.^2 > 0.65^2)= 0;

    case 1
        % Half circle membrane.
        R = 0.7;
        vD = half_circle_draw(sz, R) * velocity_scale;
        
    case 2
        % Membrane shape.
        vD = membren(sz).* velocity_scale;
        
    case 3
        % Advanced membrane (similar to real one).
        vD = advanced_membrane(sz) * velocity_scale;

    case 4
        % Advanced membrane (Different version).
        vD = advanced_membrane2(sz) * velocity_scale;

    case 5 
        % Advanced membrane with control of:
        angs = 5; % N. of angels.
        rads = 4; % N. of radiuses.
        
        vD = advanced_membrane_set(sz, angs, rads) * velocity_scale;
        
    case 6
        % Parabul membrane.
         vD = (xx.^2 + yy.^2) * velocity_scale;
         vD(xx.^2 + yy.^2 > 0.65^2)= 0;
    
    case 7
        vD = chrs_pezo * velocity_scale;
   otherwise
        disp('Invalid option for "membrane" variable');
end

% Doppler frequencies.
fD = 2.*vD./lambda; 
pcolor(xs, ys, fD);
title('Spatially-varying Velocity Pattern', 'FontSize', 13);
colorbar('eastoutside')

% PHYSICAL PARAMETERS (Raman-Nath Scattering (lamda*L / d^2 << 1): 

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
I1 = (abs(besselj(0, rho) + 1i.*besselj(1, rho).*exp(-1i.*psi)).^2); % *Ir
I0 = (abs(besselj(-1, rho) + 1i.*besselj(0, rho).*exp(-1i.*psi)).^2); % Ir

photons_after_liquid = sum(I1, 'all') + sum(I0, 'all');

balance = I1 - I0;
figure, image(balance,'CDatamapping','scaled'),set(gca,'FontSize',20);
title('Original balanced signal');
colorbar('eastoutside');

% We start with one white balance measurement.
white_bal_measre = sum(balance, 'all');

permutation1 = randperm(N, M);
permutation0 = randperm(N, M);

% Two different masks wavelets:

% Positive Dragon.
masks1 = double(dragon_masks(N, ...
       M, permutation1)');

% Positive Hadamard.
masks0 = double(hadamard_masks(N, ...
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
total_projection = y + b - white_bal_measre;
add_noise = 1;

switch add_noise
        
    case 1
        noise = std(total_projection) * (rand(M,1) - 0.5).* 1e-3;
        sigma = std(noise);
        delta = sqrt(M + 2*sqrt(2*M)) * sigma;
        
    case 0
        
        noise = 0;
end

total_projection = total_projection + noise;
% Final measurement matrix.
masks = masks1 + masks0 - 1;
masksT = masks';

% Final measurement matrix pointer function.
A = @(x)masks * x;
At = @(x)masksT * x;

AA = A; AAt = At; bb = total_projection;

delta = 0;

clear optsNESTA
optsNESTA.TypeMin = 'tv';
optsNESTA.Verbose = 10;
optsNESTA.TolVar =1e-3;
optsNESTA.stopTest = 1;
optsNESTA.xplug = masksT * total_projection; 

Type = 2;

switch Type
% four options:
% option 1: 
    case 1
 % orthogonalize A with a QR.  This works for noiseless data
        % It's meant for the case where A is a matrix, not a function
        [Q,R] = qr(masks',0);  % A = (Q*R)'= R'*Q', R is triangular
        bb = (R')\total_projection;
        AA = @(x) (R')\A(x);  % this is now Q'
        AAt= @(x) At(R\x);    % this is now Q

% option 2:
    case 2
% find a cholesky decomposition of A. Noiseless data only
        delta = 0;
        R = chol( masks*masks');
        optsNESTA.AAtinv = @(x) R\( R'\x );

% option 3:
    case 3
% Find the SVD of A.  Noisy data are OK.
            delta = 1e-5;
            [U,S,V] = svd(masks,'econ');
            optsNESTA.USV.U = U;
            optsNESTA.USV.S = S;
            optsNESTA.USV.V = V;

% option 4:
    case 4
    % Use a Krylov Subspace method to solve for inv(AA')
    % Noiseless data only.  Here, we'll use Conjugate-Gradients
    % (call CGwrapper, which calls MATLAB's "pcg" function)
    delta = 0;
    A_function = @(x) A(At(x));
    cg_tol = 1e-6; cg_maxit = 40;
    CGwrapper(); % (first, zero-out the CGwrapper counters)
    optsNESTA.AAtinv = @(total_projection) CGwrapper(A_function,b,cg_tol,cg_maxit);
    
otherwise
        disp('Invalid option for "Type" variable');
end 

muf = 1e-6; % 1e-6 seems to wor the best
[reconpic , ~] = NESTA(AA, AAt,...
   total_projection, muf, 0, optsNESTA);

diff_rec = reshape(reconpic, [sz, sz])';

figure, image(diff_rec,'CDatamapping','scaled'),set(gca,'FontSize',20);
title('Reconstructed balance NESTA');
colorbar('eastoutside');

% ERROR CALCULATIONS:

if signal_type == 1 
    % Gaussian signal with standard deviation w.
    center = xx.^2 + yy.^2 < w;
    
else 
    % Circular signal with 0.4 radius.
    center = xx.^2 + yy.^2 < sig_rad; 
    
end

% Number of elements to be taken.
flat_cen = reshape(center.',1,[])';
num = size(flat_cen(flat_cen ~= 0));
num = num(1);

rel_diff_rec = diff_rec .* center;
rel_balance = balance .* center;


% 1. Sensativity error.
all = sum(abs(balance - diff_rec), 'all');
relative_sensativity_err = all ./ sz^2;


% 1.1 Sensativity error for the relevant part.
rel_all = sum(abs(rel_balance - rel_diff_rec), 'all');
relevant_relative_sensativity_err = rel_all ./ num;

