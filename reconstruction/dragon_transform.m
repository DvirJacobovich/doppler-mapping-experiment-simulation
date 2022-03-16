function y = dragon_transform( x, picks, N, mode )
%m_transform Fast dragon noiselet transform or it's traspose. 
switch mode
    case 1
        y = A_noiselet(x,picks);
    case 2
        y = At_noiselet(x,picks, N);
    otherwise
        error('Unknown mode passed to f_handleA!');
end

end

