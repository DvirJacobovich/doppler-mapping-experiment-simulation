function y = hadamard_transform( x, picks, N, mode )
%m_transform Fast dragon noiselet transform or it's traspose. 
switch mode
    case 1
        n = length(x);
        w = hadamardc(x)/sqrt(n);
        y = w(picks);
    case 2
        vn = zeros(N,1);
        vn(picks) = x;
        y = hadamardc(vn/sqrt(N));
    otherwise
        error('Unknown mode passed to hadamard_transform!');
end

end

