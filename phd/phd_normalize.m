%normalizes phase distributions such that the sum is equal to 1. First and
%last histogram column are half width only. Adds phase distribution for 
%N=0.
function phd_normalize(input, output)

load(input, 'N', 'phis', 'profiles');

N = [0; N];

save(output, 'N', 'phis', 'profiles');

% for each profile found in 'input':
for p = 1:size(profiles, 1)
    
    temp = load(input, profiles(p, :));
    
    % add phase histogram for N=0:
    phd = [zeros(1, 1024); temp.(profiles(p, :))];
    
    % if current profile is flow compensated:
    if (profiles(p,1) == 'f')
        % for N=0 we have phi=0
        phd(1, 1) = 1;
    else
        % for N=0 equal distribution of phases,
        % first & last column basically have half width only.
        phd(1, :) = 1;
        phd(1, 1) = 0.5;
        phd(1, numel(phis)) = 0.5;
    end
    
    % normalization
    for i = 1:numel(N)
        phd(i, :) = phd(i, :) / sum(phd(i, :));
    end
    
    eval([profiles(p, :) ' = phd;']);
    
    save(output, profiles(p, :), '-append');
    
end

end