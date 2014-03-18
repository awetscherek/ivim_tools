%integrated phase distribution for basic flow compensated gradient profile 
%(FEER waveform) and N=1.
function retval = phd_int_fc0_N1(x)

    retval = 0.5 + (2 * phs(x) - 1) .* phd_int_fc0_N1_plus(abs(x));
end

%analytic formula for integrated phase distribution for basic flow
%compensated gradient profile (FEER waveform) and N=1 (only pos. phases).
function retval = phd_int_fc0_N1_plus(x)

    xmax = sqrt(3)/2;
    z = 1/4;

    retval = z/xmax * phs(xmax - x) .* phs(x) ... only for 0 < x < xmax well defined
        .* (phs(xmax/2 - x) .* (-4*x + 8*sqrt(2)/3*sqrt(xmax.*x) + 4/3*x.^2/xmax + 2*x*sqrt(2)*atanh(1/sqrt(2)) ...
        -x.^2/xmax/sqrt(2).*atanh(1/sqrt(2))-x.^2/xmax) ...
        + (1 - phs(xmax/2 - x)) ...
          .* (-2*sqrt(2)*sqrt(1-x/xmax)*xmax+2*sqrt(2)*x.*atanh(sqrt(1-max(0.00001,x/xmax))) ...
        +sqrt(2)/6*x.*sqrt(1-x/xmax)+sqrt(2)/3*xmax*sqrt(1-x/xmax)-x.^2/sqrt(2)/xmax.*atanh(sqrt(1-max(0.00001,x/xmax))) ...
        +sqrt(2)*2/3*xmax*sqrt(1- x/xmax).*(1-x/xmax) ...
        +2*xmax)) + (1 - phs(xmax - x)) .* 0.5;
end

