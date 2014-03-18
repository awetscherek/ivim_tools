%analytic expression for phase distribution in the limit of high N.
function retval = phd_lim_Nhigh(x, N)

    retval = sqrt(3*N/2/pi) * exp(-x.^2 *3*N/2);
end

