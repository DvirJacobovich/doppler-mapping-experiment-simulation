function[sens_err] = sensativity_err(balance, noiseles_Ip, sz)
    % ERROR CALCULATIONS:

%     if signal_type == 1 
%         % Gaussian signal with standard deviation w.
%         center = xx.^2 + yy.^2 < w;
% 
%     else 
%         % Circular signal with 0.4 radius.
%         center = xx.^2 + yy.^2 < sig_rad; 
% 
%     end

    % Number of elements to be taken.
%     flat_cen = reshape(center.',1,[])';
%     num = size(flat_cen(flat_cen ~= 0));
%     num = num(1);
% 
%     rel_Ip = noiseles_Ip .* center;
%     rel_balance = balance .* center;


    % 1. Sensativity error.
    all = sum(abs(balance - noiseles_Ip), 'all');
    sens_err = all ./ sz^2;


%     % 1.1 Sensativity error for the relevant part.
%     rel_all = sum(abs(rel_balance - rel_Ip), 'all');
%     relevant_relative_sensativity_err = rel_all ./ num;
    
end