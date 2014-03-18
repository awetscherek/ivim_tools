
function [tdata, tbetta, tb, tT] = trace(data, betta, b, T)
% prepares data for fitting (averaging over directions and averaging 
% b0 images within a single breath hold)

height = size(data, 2);
width = size(data, 3);

% count breath holds / acquisition series assuming that each series starts 
% and ends with an unweighted (b0) acquisition.
bh = numel(find([b;0] == [0;b] & [b;0] == 0)) - 1;
bh_size = numel(b) / bh;

% average b0 images for each breath hold and determine number of different
% diffusion weightings:
weightings = numel(find([b;0] ~= [0;b])) - numel(find(b == 0)) + bh;

% initialize output variables 
tdata = zeros(weightings, size(data, 2), size(data, 3));
tbetta = zeros(weightings, 1);
tT = tbetta;
tb = tbetta;

w_count = 0;
last = 0;

for i = 1:bh
    
    b_sub = b((i-1) * bh_size + (1:bh_size));
    
    tdata(w_count + i, :) = max(0.33, mean(data(find(b_sub == 0) + (i-1) * bh_size, :), 1));
    tT(w_count + i) = T((i-1) * bh_size + 1);
    tb(w_count + i) = 0;
    tbetta(w_count + i) = 2;
    
    first = (i-1) * bh_size + 1;
    last = i * bh_size;
    
    while (first < last)
        while (b(first) == 0)
            first = first + 1;
        end
        next = first;
        while ((next < last) && (b(next + 1) == b(first)))
            next = next + 1;
        end
        w_count = w_count + 1;
        tdata(w_count + i, :) = nthroot(prod(abs(data(first:next, :)), 1), next-first+1);
        tT(w_count + i) = T(first);
        tb(w_count + i) = b(first);
        tbetta(w_count + i) = betta(first);
        first = next + 1;
    end
end

end