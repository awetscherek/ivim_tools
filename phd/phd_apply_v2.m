%input/output phd file, N2: values for N, for which phd is calculated, 
%vsteps: discretization steps for velocity distribution.
function phd_apply_v2(input, output, N2, vsteps)

load(input, 'N', 'phis', 'profiles');

phis_old = phis;

N_old = N;
N = N2;

%maximum obtainable phase is 2x the max phase of the input distribution
phis = [phis_old phis_old(2:length(phis_old)) + max(phis_old)];

save(output, 'N', 'phis', 'profiles');

%weighting factors:
fr = ones(1, vsteps) / (vsteps - 1);
fr(1) = fr(1) / 2;
fr(vsteps) = fr(vsteps) / 2;

%relative velocities:
vr = (0:(vsteps-1)) / (vsteps-1) * 2;

%for each relative velocity find indices in original phase distribution...
iceil = max(0, ceil(log(N * vr/N_old(2)) / log(max(N_old)/N_old(2)) * (numel(N_old)-2))) + 2;
ifloor = max(-1, floor(log(N * vr/N_old(2)) / log(max(N_old)/N_old(2)) * (numel(N_old)-2))) + 2;

for p = 1:size(profiles, 1)
    
    temp = load(input, profiles(p, :));
    
    % add phase histogram for N=0:
    phd = temp.(profiles(p, :));
    
    phd_v2 = zeros(length(N), length(phis));

    for m = 1:vsteps
    
        waitbar(((m-1) / vsteps + p - 1) / size(profiles, 1));
        
        ind_high = iceil(:, m) > numel(N_old);
        ind_low = iceil(:, m) <= numel(N_old);
    
        if (vr(m) > 0)
            phd_v2(ind_high, :) = phd_v2(ind_high, :) + ...
                2 * fr(m) * phd_rebin(phis_old * vr(m), phd_hist_lim_Nhigh(phis_old, phis_old, N(ind_high) * vr(m)), phis);
            
            c = (N(ind_low) * vr(m) - N_old(ifloor(ind_low, m))) ./ (N_old(iceil(ind_low, m)) - N_old(ifloor(ind_low, m))) * ones(size(phis));
    
            phd_v2(ind_low, :) = phd_v2(ind_low, :) + ...
                fr(m) * ((1-c) .* phd_rebin(phis_old * vr(m), phd(ifloor(ind_low, m), :), phis) ...
                    + c .* phd_rebin(phis_old * vr(m), phd(iceil(ind_low, m), :), phis));
        else
           
            phd_v2(:, 1) = phd_v2(:, 1) + fr(m);
            
        end
      
    end
    
    eval([profiles(p, :) ' = phd_v2;']);
    
    save(output, profiles(p, :), '-append');

end