%replaces phase distribution for N=1 with the analytic solution for the
%gradient profiles.
function phd_exact_N1(input, output)

load(input, 'N', 'phis', 'profiles');

if (strcmp(input, output))
    save(output, 'N', 'phis', 'profiles', '-append');
else
    save(output, 'N', 'phis', 'profiles');   
end

if (~isempty(find(N == 1, 1)))
    
    ind = find(N == 1, 1);

    % for each profile found in 'input':
    for p = 1:size(profiles, 1)
    
        temp = load(input, profiles(p, :));
       
        % add phase histogram for N=0:
        phd = temp.(profiles(p, :));
    
        eval(['phd(ind, :) = phd_hist_' profiles(p, :) '_N1(phis, phis);']);
        phd(ind, :) = phd(ind, :) / sum(phd(ind, :));
        
        eval([profiles(p, :) ' = phd;']);
    
        save(output, profiles(p, :), '-append');

    end  
    
end

end