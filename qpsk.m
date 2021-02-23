function d = qpsk(b)
% d = qpsk(b)
%
% Map the bits to be transmitted into QPSK symbols using Gray coding. The
% resulting QPSK symbol is complex-valued, where one of the two bits in each
% QPSK symbol affects the real part (I channel) of the symbol and the other
% bit the imaginary part (Q channel). Each part is subsequently PAM
% modulated to form the complex-valued QPSK symbol. The energy per QPSK
% symbol is normalized to unity.
%
% The mapping resulting from the two PAM branches are:
%
% complex part (Q channel)
%         ^
%         |
%  10 x   |   x 00   (odd bit, even bit)
%         |
%  -------+------->  real part (I channel)
%         |
%  11 x   |   x 01
%         |
%
%
%
% Input:
%   b = bits {0, 1} to be mapped into QPSK symbols
%
% Output:
%   d = complex-valued QPSK symbols

%Just for simulation and cross checking, comment before running the final
%code
%nr_data_bits=1000;
%b_data = (randn(1, nr_data_bits) > .5);
%b = [b_data];

d = zeros(1, length(b)/2);

for n = 1:length(b)/2
    odd = b(2*n);
    even = b(2*n-1);
    if (even == 0) && (odd == 0)
        d(n) = exp(1i*pi/4);   % This denotes 45'
    end
    if (even == 1) && (odd == 0)
        d(n) = exp(1i*3*pi/4); % This denotes 135'
    end
    if (even == 1) && (odd == 1)
        d(n) = exp(1i*5*pi/4); % This denotes 225'
    end
    if (even == 0) && ( odd== 1)
        d(n) = exp(1i*7*pi/4); % This denotes 315'
    end
end

figure;
plot(d,'x');
xlabel('Real'); ylabel('Imaginary');
title('QPSK constellation');

