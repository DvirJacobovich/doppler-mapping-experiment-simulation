
clear
close all force

addpath('C:\Users\Dvir\Desktop\Doppler Frequencies Mapping')
addpath('C:\Users\Dvir\Desktop\Doppler Frequencies Mapping\Compressive_Sensing_Tools\l1magic\Optimization')
addpath 'C:\Users\Dvir\Desktop\Doppler Frequencies Mapping\Usefull stuff'
addpath 'C:\Users\Dvir\Desktop\DMD Settings and Control\Compressive_Sensing_Tools\NESTA_v1.1\Analysis'


addpath utils/
addpath utils/qr/
addpath utils/utils_Wavelet

% No. of reweighting iterations (Optional).
rwt = 5;       
pow = 6;
sz = 2^pow;

N = sz^2; 

% No. of measurements
M = round((N)/10); 

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

% Liquid crystal relaxation time.
tau = 72*10^-3; %[s].

% Gauss signal sqrt(std).
w = 0.4;

% Uniformed signal radius.
sig_rad = 0.4;

% Signal amplitude [mW / cm^2].
Amp = 5; 

t = 1;
train_data_struct = struct('sample_num', 0, 'balance', zeros(sz, sz), 'noisy_rec', zeros(sz, sz),...
        'noisy_err', 0, 'noiseles_rec', zeros(sz, sz), 'noiseles_err', 0);

i = 1; trian_data_len = 2;

while i <= trian_data_len
    
    try    
        signal_type = randi(2);
        membrane_type = randi([2 5]);

        switch signal_type
            case 1
                % Gaussian signal with incertainty of sqrt(w), [sqrt(mW / cm^2)].
                r = Amp.*exp(-((xx.^2) + (yy.^2)) ./ w).*exp(-1i.*pi*t); 

            case 2
                % Circular signal radius sig_rad, [sqrt(mW / cm^2)].
                r = Amp.*circle_draw(sz, sig_rad); 
        end
        
        % Signal intensity in [mW / cm^2].
        Ir = abs(r).^2; 

        % Beam wavelenght.
        lambda = 532e-9; % [m].

        % Piazo engeen doppler frequencies membren.
        velocity_scale = 1e-8; % [m / s].

        switch membrane_type
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
                vD = membren(sz) * velocity_scale;

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

           otherwise
                disp('Invalid option for "membrane" variable');

        end
        pcolor(xs, ys, vD);
        title('Spatially-varying Velocity Pattern', 'FontSize', 13);
        colorbar('eastoutside')

        % Doppler frequencies.
        fD = 2.*vD./lambda; 

        rho = (2*k0*n2.*d ./ (((1 + ld^2*kg^2)^2 + (fD.*tau)^2).^0.5)).*(Ir.*Ir).^0.5;
        psi = atan(fD.*tau ./ (1 + ld^2*kg^2)^2);

        % After liquid crystal 
        I1 = (abs(besselj(0, rho) + 1i.*besselj(1, rho).*exp(-1i.*psi)).^2).*Ir; %/Ir
        I0 = (abs(besselj(-1, rho) + 1i.*besselj(0, rho).*exp(-1i.*psi)).^2).*Ir;% /Ir

        balance = I1 - I0;
        figure, image(balance,'CDatamapping','scaled'),set(gca,'FontSize',20);
        title('Original balanced signal');
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
        masks0 = double(hadamard_masks(N, ...
               M, permutation0)');


        b1 = masks1 * reshape(I1.',1,[])';
        b0 = masks0 * reshape(I0.',1,[])';

        % First balanced projections.
        b = b1 - b0; 

        y1 = masks0 * reshape(I1.',1,[])';
        y0 = masks1 * reshape(I0.',1,[])';

        % Second balanced projections.
        y = y1 - y0; 

        % Subtracting the White Balanced measurment.
        total_projection = y + b - white_bal_measre;
        
        % Noise calculation.
        noise = std(total_projection) * (rand(M,1) - 0.5).* 1e-2;
        sigma = std(noise); delta = sqrt(M + 2*sqrt(2*M)) * sigma;
        noisy_tot_proj = total_projection + noise;  
        muf = 1e-10;

        % Final measurement matrix.
        A = masks1 + masks0 - 1;
        At = A';

        % Final measurement matrix pointer function.
        A_f = @(x)A * x;
        At_f = @(x)At * x;

        % Initial gues.
        x0 = At * total_projection; 

        tvI = sum(sum(sqrt([diff(balance,1,2) zeros(sz,1)].^2 + ...
                                   [diff(balance,1,1); zeros(1,sz)].^2 )));
                               
        fprintf('Original TV = %.3f\n', tvI);
        time0 = clock;

        % OPTIMIZATION PARAMETERS:

        % Duality gap barrier - Algorithm terminates when the duality gap <= lbtol.
        lbtol = 1e-3; % has been changed frm 1e-3
        % Tolerance for SYMMLQ.
        slqtol = 1e-8;
        % Factor by which to increase the barrier constant at each iteration.
        mu = 10;
        % slqmaxiter - Maximum number of iterations for SYMMLQ;
        slqmaxiter = 200;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%% noiseles part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        noiseles_rec =  tveq_logbarrier(x0, A_f, At_f, total_projection, lbtol, mu,...
            slqtol, slqmaxiter); 

        noiseles_Ip = reshape(noiseles_rec, sz, sz)';
        figure, image(noiseles_Ip,'CDatamapping','scaled'),set(gca,'FontSize',20);
        title('Noiseles balance rec. Per 100%');
        colorbar('eastoutside')
        fprintf('Total elapsed time noiseles = %f secs\n\n', etime(clock,time0));

        noiseles_sens_err = sensativity_err(balance, noiseles_Ip, sz);

        %%%%%%%%%%%%%%%%%%%%%%%%% noisy part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        noisy_rec =  tveq_logbarrier(x0, A_f, At_f, noisy_tot_proj, lbtol, mu,...
            slqtol, slqmaxiter); 

        noisy_Ip = reshape(noisy_rec, sz, sz)';
        figure, image(noisy_Ip,'CDatamapping','scaled'),set(gca,'FontSize',20);
        title('Noisy balance rec. Per 100%');
        colorbar('eastoutside')
        fprintf('Total elapsed time noisy = %f secs\n\n', etime(clock,time0));

        noisy_sens_err = sensativity_err(balance, noisy_Ip, sz);

        s_i = struct('sample_num', i, 'balance', balance, 'noisy_rec', noisy_Ip,...
            'noisy_err', noisy_sens_err, 'noiseles_rec', noiseles_Ip, 'noiseles_err', ...
            noiseles_sens_err);

        train_data_struct(end + 1 : end + length(s_i)) = s_i;        
        struct_name = sprintf('train_data_struct_batch_number_%d.mat', i);
        
        save(struct_name, 'train_data_struct');
        
        if mod(i, 5) == 0
            struct_name = sprintf('train_data_struct_batch_number_%d.mat', i / 5);
            train_data_struct = struct('sample_num', 0, 'balance', zeros(sz, sz),...
                'noise', noise, 'noisy_rec', zeros(sz, sz), 'noisy_err', 0, ...
                'noiseles_rec', zeros(sz, sz), 'noiseles_err', 0);

            save(struct_name, 'train_data_struct');
        end

    catch
        
        fprintf('Couldnt reconstruct try again with different measurement matrix');
        continue;
        
    end
    
    i = i + 1;
    
end


                                                                                                                         
