%integrated phase distribution for limit of high N.
function retval = phd_int_lim_Nhigh(x, N)

    retval = (1 + erf(sqrt(3 * N / 2) * x)) / 2;

end

