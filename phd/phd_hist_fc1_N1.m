%calculates histogram for phase distribution for the symmetric flow 
%compensated gradient profile and N=1.
function retval = phd_hist_fc1_N1(x, phis)
    retval = phd_int_fc1_N1(min(max(phis), abs(x) + phis(2) / 2)) ...
        - phd_int_fc1_N1(max(min(phis), abs(x) - phis(2) / 2));
end

