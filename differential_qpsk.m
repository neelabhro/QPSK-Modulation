function d = differential_qpsk(b)
d = qpsk(b);
n = length(d);
for sym = 2:n
    d(sym) = d(sym)*d(sym-1);
end