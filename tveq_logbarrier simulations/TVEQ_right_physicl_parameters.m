    % close all force
    
   addpath 'C:\Users\Dvir\Desktop\Doppler Frequencies Mapping'
   addpath 'C:\Users\Dvir\Desktop\Doppler Frequencies Mapping\Compressive_Sensing_Tools\l1magic\Optimization'
    addpath '.../'
      
    pow = 7;
    sz = 2^pow;

    N = sz^2; 

    % No. of measurements
    M = round(N/8); 

    xs = linspace(-1, 1, sz);
    ys = xs;

    [xx, yy] = meshgrid(xs, ys);

    % PHYSICAL PARAMETERS (Raman-Nath Scattering approx (lamda*L / d^2 << 1): 

    % Beams wave lenghts, in [m].
    lambda = 532e-9; 

    % Magnitude of the laser wavevector.
    k0 = 2*pi / lambda; % in [1/m]

    % Kerr like Refractive index.
    n2 = -1.8*10^-4; % [m^2 / mW].

    % Fringe spacing.
    d = 50*10^-6; % [m].

    % Interaction region thickness.
    L = 9*10^-6; % [m]. 

    ld = 12*10^-6; % [m].

    % Grating index wavevector.
    kg = 2*pi/(110*10^-6); % [1 / m].

    % Liquid crystal relaxation time under these parameters.
    tau = 72*10^-3; %[s].

    % Gauss signal sqrt(std).
    w = 0.4;

    % Uniformed signal radius.
    sig_rad = 0.4;

    t = 1;
  
    % In the future the fD piezo will take place of the balance.
    % train_data_struct = struct('sample_num', 0, 'balance', zeros(sz, sz),...
    %       'meaures', zeros(sqrt(M), sqrt(M)), 'piezo', zeros(sz));
    
    % Hz frequencies linear range.
    Hz_freqs = linspace(0.3, 2, 10);
    v_scale = lambda .* Hz_freqs;
    % Chosen freq.
    Hz = Hz_freqs(randi(10));
    
    % Rad/sec freq.
    omega = 2*pi*Hz;
    
    % Piezo velocity.
    v_scale = lambda * Hz;
    % lambda * Hz = v = A*omega*cos(omega*t) => in t = 0 we have A*omega =
    % Hz * lambda => A = (Hz * lambda) / omega.
    
    % Signal amplitude [V / cm^2].
    Amp = sqrt((Hz * lambda) / omega); % I added sqrt since we want the 
    i = 1; % trian_data_len = 1;

        signal_type = 1;

        switch signal_type
            case 1
                % Gaussian signal with incertainty of sqrt(w).
                r = Amp*exp(-((xx.^2) + (yy.^2)) ./ w).*exp(-1i.*pi*t);

            case 2
                % Circular signal radius sig_rad.
                r = Amp*circle_draw(sz, sig_rad);
        end

    Ir = abs(r).^2; 

    % Refferenced signal (has no usge now).
    s = r;
    Is = abs(s).^2; 

    % Beams wave lenghts, in [m].
    lambda = 532e-9;

    % Piazo engeen doppler frequencies membren.
    Hz_rage = linspace(-0.5, 0.5, 10);
    
    % Doppler frequencies.
    piezo_freqs = last_cris_piexo(sz);
    figure, pcolor(xs, ys, piezo_freqs);
    title(sprintf('Spacialy varying freqs of %d Piezos [-0.8Hz 0,8Hz] raged', 40));
    colorbar('eastoutside')
    
    % Piazo engeen velocities.
    piezo_vels = piezo_freqs .* lambda;
    figure, pcolor(xs, ys, piezo_vels);
    title(sprintf('Spacialy varying velocities of %d Piezos [m/s]', 40));
    colorbar('eastoutside')

    % PHYSICAL PARAMETERS (Raman-Nath Scattering (lamda*L / d^2 << 1): 

    rho = (2*k0*n2.*d ./ (((1 + ld^2*kg^2)^2 + (piezo_freqs.*tau)^2).^0.5)).*(Is.*Ir).^0.5;
    psi = atan(piezo_freqs.*tau ./ (1 + ld^2*kg^2)^2);

    % After liquid crystal 
    I1 = (abs(besselj(0, rho) + 1i.*besselj(1, rho).*exp(-1i.*psi)).^2); % *Ir
    I0 = (abs(besselj(-1, rho) + 1i.*besselj(0, rho).*exp(-1i.*psi)).^2); % Ir

    balance = I1 - I0;
    figure, image(balance,'CDatamapping','scaled'),set(gca,'FontSize',20);
    title('Balance signal');
    colorbar('eastoutside')
    
    % We start with one white balance measurement.
    white_bal_measre = sum(balance, 'all');

    permutation1 = randperm(N, M);
    permutation0 = randperm(N, M);

    % Two different masks wavelets:

    % Positive Dragon.
    masks1 = double(dragon_masks(N, ...
           M, permutation1)');

    % Positive Hadamard.
    masks0 = double(dragon_masks(N, ...
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
    
    add_noise = 0;
    switch add_noise 
        case 1 
            noise = std(total_projection) * (rand(M,1) - 0.5).* 1e-3;
            sigma = std(noise);
            delta = sqrt(M + 2*sqrt(2*M)) * sigma;
            tot_proj = total_projection + noise;
    
        case 0
            tot_proj = total_projection;
    end
   

    
   % Final measurement matrix.
A = masks1 + masks0 - 1;
At = A';

% Final measurement matrix pointer function.
A_f = @(x)A * x;
At_f = @(x)At * x;

% Initial gues.
x0 = At * total_projection; 

tvI = sum(sum(sqrt([diff(balance,1,2) zeros(sz,1)].^2 + [diff(balance,1,1); zeros(1,sz)].^2 )));
fprintf('Original TV = %.3f\n', tvI);

time0 = clock;

% OPTIMIZATION PARAMETERS:

% Duality gap barrier - Algorithm terminates when the duality gap <= lbtol.
lbtol = 1e-11; % has been changed frm 1e-3
% Tolerance for SYMMLQ.
slqtol = 1e-8;
% Factor by which to increase the barrier constant at each iteration.
mu = 10;
% slqmaxiter - Maximum number of iterations for SYMMLQ;
slqmaxiter = 200;

xp =  tveq_logbarrier(x0, A_f, At_f, total_projection, lbtol, mu, slqtol, slqmaxiter); 

Ip = reshape(xp, sz, sz)';
figure, image(Ip,'CDatamapping','scaled'),set(gca,'FontSize',20);
title('Reconstruction balanced Per 12.5%');
colorbar('eastoutside')

    
