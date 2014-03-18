%weighting function for bipolar gradient profile
function retval = c_bip(r, N)
    
    lb = min(1, max(0, (r - 1)./N));
    ub = min(1, r./N);
    
    retval = 2 * sqrt(3) * (-ub.^2/2 + max(0, ub - 0.5) .* (ub - 0.5) ...
        + lb.^2/2 - max(0, lb - 0.5) .* (lb - 0.5));
    
end