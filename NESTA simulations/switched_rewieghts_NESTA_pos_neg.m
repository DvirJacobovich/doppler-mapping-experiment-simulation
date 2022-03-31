clear
close all force

addpath 'C:\Users\Dvir\Desktop\DMD Settings and Control\Compressive_Sensing_Tools\NESTA_v1.1\Analysis'


addpath utils/
addpath utils/qr/
addpath utils/utils_Wavelet

% No. of reweighting iterations (Optional).
pow = 6;
sz = 2^pow;

N = sz^2; 

% No. of measurements
M = round(N); 

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

% The smaller the velocity scale is, the smaller the balanced scale.
velocity_scale = 1e-9; 
R = 0.7;
vD = half_circle_draw(sz, R) * velocity_scale;

pcolor(xs, ys, vD);
title('Spatially-varying Velocity Pattern', 'FontSize', 13);

% Beam wavelenght.
lambda = 532e-9;

% Doppler frequencies.
fD = 2.*vD./lambda; 

% PHYSICAL PARAMETERS (Raman-Nath Scattering (lamda*L / d^2 << 1): 

% Beams wave lenghts, in [m].
lambda = 532e-9; 

% Magnitude of the laser wavevector, in [1/m].
k0 = 2*pi / lambda; 

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

rho = (2*k0*n2.*d ./ (((1 + ld^2*kg^2)^2 + (fD.*tau)^2).^0.5)).*(Ir.*Ir).^0.5;
psi = atan(fD.*tau ./ (1 + ld^2*kg^2)^2);

% After liquid crystal 
I1 = (abs(besselj(0, rho) + 1i.*besselj(1, rho).*exp(-1i.*psi)).^2); % *Ir
I0 = (abs(besselj(-1, rho) + 1i.*besselj(0, rho).*exp(-1i.*psi)).^2); % Ir

% photons_after_liquid = sum(I1, 'all') + sum(I0, 'all');
% 
balance = I1 - I0;
figure, image(balance,'CDatamapping','scaled'),set(gca,'FontSize',20);
title('Original balanced signal');

% We start with one white balance measurement.
white_bal_measre = sum(balance, 'all');

not_positive_definite = 1;

while not_positive_definite
    
    try 

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
            noise_ratio = 1e-3;
            epsilon = noise_ratio .* total_projection;
            total_projection = total_projection + epsilon;

            % Not relevant since delta = 0 in the alg.
            sigma = std(epsilon);
            delta = sqrt(M + 2*sqrt(2*M))*sigma;
            muf = 1e-6 * noise_ratio;

        case 0
            delta = 0;
            muf = 1e-6;

    end

    % Final measurement matrix.
    masks = masks1 + masks0 - 1;
    masksT = masks';

    % Final measurement matrix pointer function.
    A = @(x)masks * x;
    At = @(x)masksT * x;

    AA = A; AAt = At; bb = total_projection;

    clear optsNESTA
    optsNESTA.TypeMin = 'tv';
    optsNESTA.Verbose = 10;
    optsNESTA.TolVar =1e-3;
    optsNESTA.stopTest = 2;
    optsNESTA.xplug = masksT * total_projection; 

    Type = 2;

    switch Type
    % Four options:

        case 1
        % orthogonalize A with a QR.  This works for noiseless data
        % It's meant for the case where A is a matrix, not a function

            [Q,R] = qr(masks',0);  % A = (Q*R)'= R'*Q', R is triangular
            bb = (R')\total_projection;
            AA = @(x) (R')\A(x);  % this is now Q'
            AAt= @(x) At(R\x);    % this is now Q

        case 2
        % Find a cholesky decomposition of A. Noiseless data only.
             R = chol( masks*masks');
             optsNESTA.AAtinv = @(x) R\( R'\x );

        case 3
        % Find the SVD of A.  Noisy data are OK.
            delta = 1e-3;
            [U,S,V] = svd(masks,'econ');
            optsNESTA.USV.U = U;
            optsNESTA.USV.S = S;
            optsNESTA.USV.V = V;

        case 4
        % Use a Krylov Subspace method to solve for inv(AA')
        % Noiseless data only.  Here, we'll use Conjugate-Gradients
        % (call CGwrapper, which calls MATLAB's "pcg" function)
            delta = 0;
            A_function = @(x) A(At(x));
            cg_tol = 1e-6; cg_maxit = 40;
            CGwrapper(); % (first, zero-out the CGwrapper counters)
            optsNESTA.AAtinv = @(total_projection) CGwrapper(A_function,total_projection,cg_tol,cg_maxit);

    otherwise
            disp('Invalid option for "Type" variable');
    end 
    
        not_positive_definite = 0;
        break;
        
    catch
        
        continue;
        
    end
    
end

% muf = 1e-6; % 1e-6 seems to wor the best
[reconpic , ~] = NESTA(AA, AAt,...
   total_projection, muf, 0, optsNESTA);

diff_rec = reshape(reconpic, [sz, sz])';

figure, image(diff_rec,'CDatamapping','scaled'),set(gca,'FontSize',20);
title('Reconstructed balance NESTA');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

optsNESTA.U = diag(1 ./ abs(reconpic));
optsNESTA.Ut = optsNESTA.U;
optsNESTA.normU = max( 1./(abs(reconpic)+.1) );


for i=1:10
    [reconpic , ~] = NESTA(AA, AAt,...
    bb, 0, delta, optsNESTA); 
    optsNESTA.U = diag(1 ./ reconpic);
    optsNESTA.Ut = optsNESTA.U;
    optsNESTA.normU = max( 1./(abs(reconpic)+.1) );
    optsNESTA.xplug = reconpic;  % use old solution as starting value

end
fprintf('l2 error: %.2e\n',norm(reshape(reconpic, [sz, sz]) - balance )/norm(balance));
% reconpic = generate_haar(sz^2)' * reconpic; % for sparse rep rec.
diff_rec = reshape(reconpic, [sz, sz])';

% norm to [1,0].

figure, image(diff_rec,'CDatamapping','scaled'),set(gca,'FontSize',20);
title('Reconstructed diff');

% Error calculations.

% 1. General relative (per pixel) error.
D = abs(diff_rec - balance).^2;
relative_MSE_err = sum(D(:)) / numel(diff_rec);


% 2. Only relevant area error.
if signal == 1 
    % Gaussian signal with standard deviation w.
    center = xx.^2 + yy.^2 < w;
    
else 
    % Circular signal with 0.4 radius.
    center = xx.^2 + yy.^2 < 0.4; 
end

rel_Ip = diff_rec .* center;
rel_balance = balance .* center;

% Number of elements to be taken.
flat_cen = reshape(center.',1,[])';
num = size(flat_cen(flat_cen ~= 0));
num = num(1);

D_rel = abs(rel_Ip - rel_balance).^2;
relative_rel = sum(D_rel(:)) / num;


% in case where fD*tau <<(1 + (ld*kg)^2) we get this approximation:

function[norm_mat] = one_zore_norm(mat)
    min_mat = min(min(mat));
    max_mat = max(max(mat));
    norm_mat = (mat - min_mat) ./ (max_mat - min_mat);
    
end
