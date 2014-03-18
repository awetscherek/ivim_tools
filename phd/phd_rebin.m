%rebinning of histogram. Important: both input and output histogram have to
%cover the same bin range and have to have constant bin widths.
function hist2 = phd_rebin(phis1, hist1, phis2)

M = zeros(length(phis2), length(phis1));

if (phis1(2) > phis2(2))
   
    left = ceil((phis2 - phis2(2)/2 + phis1(2)/2) / phis1(2));
    right = ceil((phis2 + phis2(2)/2 + phis1(2)/2) / phis1(2));
    
    right = right(right <= numel(phis1));
    left = left(left <= numel(phis1));
    
    M((1:numel(right)) + left(1:numel(right)) * numel(phis2) - numel(phis2)) = ...
        M((1:numel(right)) + left(1:numel(right)) * numel(phis2) - numel(phis2)) + ...
            (phis1(left(1:numel(right))) + phis1(2)/2 - phis2(1:numel(right)) + phis2(2)/2) / phis1(2);
    M((1:numel(right)) + right * numel(phis2) - numel(phis2)) = ...
        M((1:numel(right)) + right * numel(phis2) - numel(phis2)) + ...
            (-phis1(left(1:numel(right))) - phis1(2)/2 + phis2(1:numel(right)) + phis2(2)/2) / phis1(2);
        
    M(2) = 2 * M(2);
        
    if (numel(left) > numel(right))
        M(numel(left) + left(numel(left)) * numel(phis2) - numel(phis2)) = ...
            M(numel(left) + left(numel(left)) * numel(phis2) - numel(phis2)) + ...
                (phis1(left(numel(left))) + phis1(2)/2 - phis2(numel(left)) + phis2(2)/2) / phis1(2);
    end
    
else
    
    left = ceil((phis1 - phis1(2)/2 + phis2(2)/2) / phis2(2));
    right = ceil((phis1 + phis1(2)/2 + phis2(2)/2) / phis2(2));
    
    M(left + (1:numel(left)) * numel(phis2) - numel(phis2)) = ...
        M(left + (1:numel(left)) * numel(phis2) - numel(phis2)) + ...
            (phis2(left) + phis2(2)/2 - phis1 + phis1(2)/2) / phis1(2);
    M(right + (1:numel(right)) * numel(phis2) - numel(phis2)) = ...
        M(right + (1:numel(right)) * numel(phis2) - numel(phis2)) + ...
            (-phis2(left) - phis2(2)/2 + phis1 + phis1(2)/2) / phis1(2);
    
end

hist2 = (M * hist1')';

end

