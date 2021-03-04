function bhat = differential_detect(r_d)
r = zeros(1, length(r_d));
r(1) = r_d(1);

n = length(r_d);
for i = 2:n
    r(i)=r_d(i)/r_d(i-1);
end
    
bhat = zeros(1,2*n);
for i= 1:n
    if (real(r(i))>0 && imag(r(i))>0)     % This denotes 00
        bhat(2*i-1) = 0;
        bhat(2*i) = 0;
    elseif (real(r(i))<0 && imag(r(i))>0) % This denotes 10
        bhat(2*i-1) = 1;
        bhat(2*i) = 0;
    elseif (real(r(i))<0  && imag(r(i))<0) % This denotes 11
        bhat(2*i-1) = 1;
        bhat(2*i) = 1;
    else                                   % This denotes 01
        bhat(2*i-1) = 0;
        bhat(2*i) = 1;
    end
end

