function phihat = phase_estimation(r, b_train)
% phihat = phase_estimation(r, b_train)
%
% Phase estimator using the training sequence. The phase estimate is
% obtained by minimizing the norm of the difference between the known
% transmitted QPSK-modulated training sequence and the received training
% part. NB! There are other ways of estimating the phase, this is just
% one example.
%
% Input:
%   r       = received baseband signal
%   b_train = the training sequence bits
%
% Output:
%   phihat     = estimated phase

phihat = 0;
qpsk_b = qpsk(b_train);
r_bits = r(1:length(qpsk_b));
e = norm(r_bits-qpsk_b);
%Minimizing the distance between the QPSK training sequence and the initial
%bits of the sync signal
for i = -3.1416:0.01:3.1416
    rprime = r_bits*exp(-1i*i);
    eprime = norm(rprime - qpsk_b);
    if(eprime < e)
        e = eprime;
        phihat = i;
    end
end
