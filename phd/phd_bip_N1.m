%analytic formula for phase distribution for bipolar gradient profile and
%N=1.
function retval = phd_bip_N1(x)
    
retval = phs(sqrt(3) / 2 - abs(x)) * sqrt(2/3).* ...
         ((1/2 + abs(x)/sqrt(3)) .* atanh(sqrt(1/2 - abs(x)/sqrt(3))) ...
         + (1/2 - abs(x)/sqrt(3)) * (atanh(1/sqrt(2)) - sqrt(2)) ...
         + sqrt(1/2 - abs(x)/sqrt(3)));
end

