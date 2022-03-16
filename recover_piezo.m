function[err, centerized_rec, centerized_orig] = recover_piezo(Ir, balance, rec_balance, fD, sz, plots)

% PHYSICAL PARAMETERS (Raman-Nath Scattering approx (lamda*L / d^2 << 1): 

% System's resolutions
xs = linspace(-1, 1, sz);
ys = xs;

[xx, yy] = meshgrid(xs, ys);

% Beams sqrt(std).
w = 0.4;

% Beams wave lenghts, in [m].
lambda = 532e-9; 

% Magnitude of the laser wavevector.
k0 = 2*pi / lambda; % in [1/m]

% Kerr like Refractive index.
n2 = -1.8*10^-4; % [m^2 / mW].

% Fringe spacing.
d = 50*10^-6; % [m].

% The linear factor.
xi = 8*pi.*besselj(0, 2*k0*d*n2.*Ir).* besselj(1, 2*k0*d*n2*Ir);
fD_rep = balance ./ (Ir .* xi);

% TODO - there is a units problem so that the inverse with different scale
% than the orig. In the meanwhile I solve it by saving the factor in div, 
% then multiplying fD_rec by div. Solve this problem.

div = fD ./ fD_rep;
div(isnan(div))=0;

fD_rec = rec_balance ./ (Ir .* xi); % This is not a good solution!
fD_rec = fD_rec .* div;


% Centerizing the orig/rec piazo by the beam std - w, since the
% reconstruction has relative big error in the edges of w (negate values).

center = xx.^2 + yy.^2 < w;

centerized_rec = fD_rec .* center;
% centerized_rec = rec_balance .* center;
centerized_orig = fD .* center;

% Showing centerizing piazos 
if plots == 1
    figure, pcolor(xs, ys, centerized_orig);
    title('Centerized Original Piezo');
    colorbar('eastoutside');

    figure, pcolor(xs, ys, centerized_rec);
    title('Centerized Recovered Piezo');
    colorbar('eastoutside');
end

% Calculating error:

% Number of elements to be taken.
flat_cen = reshape(center.',1,[])';
num_elements = size(flat_cen(flat_cen ~= 0));
num_elements = num_elements(1);

err = sum((abs(centerized_orig - centerized_rec)), 'all') ./num_elements;

end 