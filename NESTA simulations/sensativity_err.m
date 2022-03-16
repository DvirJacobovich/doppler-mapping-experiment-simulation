function[sens_err] = sensativity_err(balance, noiseles_Ip, sz)

    % 1. Sensativity error.
    all = sum(abs(balance - noiseles_Ip), 'all');
    sens_err = all ./ sz^2;


%     % 1.1 Sensativity error for the relevant part.
%     rel_all = sum(abs(rel_balance - rel_Ip), 'all');
%     relevant_relative_sensativity_err = rel_all ./ num;
    
end