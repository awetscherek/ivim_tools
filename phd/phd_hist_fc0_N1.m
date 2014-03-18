%calculates histogram for phase distribution for basic flow compensated
%gradient profile (FEER waveform) and N=1.
function retval = phd_hist_fc0_N1(x, phis)
    retval = phd_int_fc0_N1(min(max(phis), abs(x) + phis(2) / 2)) ...
        - phd_int_fc0_N1(max(min(phis), abs(x) - phis(2) / 2));
end

