%analytic formula for phase distribution for symmetric flow compensated
%gradient profile for N=1.
function retval = phd_fc1_N1(x)

    xmax = sqrt(9/4*sqrt(2) - 3);
    z = 1/2 - sqrt(2)/4;

    retval = z/xmax * phs(xmax - abs(x)) ... only for abs(x) < xmax well defined
        .* (phs(xmax/2 - abs(x)) .* (-4 + 4*sqrt(2)/3*sqrt(xmax./abs(x)) + 8/3*abs(x)/xmax + 2*sqrt(2)*atanh(1/sqrt(2)) ...
        -abs(x)/xmax*sqrt(2).*atanh(1/sqrt(2))-2*abs(x)/xmax ...
        + 2*sqrt(2)*atanh(sqrt(1-abs(x)/xmax)) -sqrt(2)*abs(x)/xmax.*atanh(sqrt(1-abs(x)/xmax)) - sqrt(2)*sqrt(1- abs(x)/xmax)) ...
        + (1 - phs(xmax/2 - abs(x))) * 2 ...
          .* (2*sqrt(2)*atanh(sqrt(1-abs(x)/xmax)) -sqrt(2)*abs(x)/xmax.*atanh(sqrt(1-abs(x)/xmax)) - sqrt(2)*sqrt(1- abs(x)/xmax)));
end

