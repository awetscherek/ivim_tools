%analytic formula for phase distribution for basic flow compensated 
%gradient profile (FEER waveform) and N=1.
function retval = phd_fc0_N1(x)

    xmax = sqrt(3)/2;
    z = 1/4;

    retval = z/xmax * phs(xmax - abs(x)) ... only for abs(x) < xmax well defined
        .* (phs(xmax/2 - abs(x)) .* (-4 + 4*sqrt(2)/3*sqrt(xmax./abs(x)) + 8/3*abs(x)/xmax + 2*sqrt(2)*atanh(1/sqrt(2)) ...
        -abs(x)/xmax*sqrt(2).*atanh(1/sqrt(2))-2*abs(x)/xmax) ...
        + (1 - phs(xmax/2 - abs(x))) ...
          .* (2*sqrt(2)*atanh(sqrt(1-abs(x)/xmax)) -sqrt(2)*abs(x)/xmax.*atanh(sqrt(1-abs(x)/xmax)) - sqrt(2)*sqrt(1- abs(x)/xmax)));
end

