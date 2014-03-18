%weighting function for simple flow compensated gradient profile
%(FEER waveform).
function retval = c_fc0(r, N)
    
    lb = min(1, max(0, (r - 1)./N));
    ub = min(1, r./N);
    
    retval = 4 * sqrt(3) * (ub.^2/2 - max(0, ub - 0.25) .* (ub - 0.25) ...
        + max(0, ub - 0.75) .* (ub - 0.75) - lb.^2/2 ...
        + max(0, lb - 0.25) .* (lb - 0.25) - max(0, lb - 0.75) .* (lb - 0.75));

end