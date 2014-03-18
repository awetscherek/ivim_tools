%calculates histogram for phase distribution in the limit of high N.
function retval = phd_hist_lim_Nhigh(x, phis, N)
    retval = phd_int_lim_Nhigh(min(max(phis), abs(x) + phis(2) / 2), N) ...
        - phd_int_lim_Nhigh(max(min(phis), abs(x) - phis(2) / 2), N);
end
