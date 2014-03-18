
function [tdata, tbetta, tb, tT] = trace14(data, betta, b, T)
% prepares data for fitting (averaging over directions and averaging
% b0 images within a single breath hold)

height = size(data, 2);
width = size(data, 3);

% only 1 b0 image for each breath hold.
bh_count = 1;
for i = 2:length(b)
    if (b(i) + b(i-1) == 0)
        bh_count = bh_count + 1;
    end
end

% average b0 images for each breath hold and determine number of different
% diffusion weightings:
weightings = 0;

    averagedb0 = zeros(bh_count, height, width);
    last = 0;
    for i = 1:bh_count
        first = last + 1;
        if (i < bh_count)
            last = first - 1 + find(b(first:(length(b)-1)) + b((first+1):length(b)) == 0, 1, 'first');
        else
            last = length(b);
        end
        temp = b(b(first:last) > 0);
        weightings = weightings + numel(find(temp(2:length(temp)) ~= temp(1:(length(temp)-1))));
        temp2 = abs(data(first:last, :, :));
        averagedb0(i, :, :) = reshape(mean(temp2(b(first:last) == 0, :, :), 1), 1, height, width);
    end


tdata = zeros(weightings, size(data, 2), size(data, 3));
tbetta = zeros(weightings, 1);
tT = tbetta;
tb = tbetta;

w_count = 0;

last = 0;
for i = 1:bh_count
    tempb0 = reshape(averagedb0(i, :, :), height, width);
    first = last + 1;
    last = length(b); 

        if (i < bh_count)
            last = first - 1 + find(b(first:(length(b)-1)) + b((first+1):length(b)) == 0, 1, 'first');
        end

    tdata(w_count + i, :, :) = max(0.33, tempb0);
    tT(w_count + i) = T(first);
    tb(w_count + i) = 0;
    tbetta(w_count + i) = 2;
    
    while (first < last)
        while (b(first) == 0)
            first = first + 1;
        end
        next = first;
        while ((next < last) && (b(next + 1) == b(first)))
            next = next + 1;
        end
        tempb = reshape(sum(abs(data(first:(next-2), :, :)), 1) / (next-first-1), height, width);
        w_count = w_count + 1;
        tdata(w_count + i, :, :) = tempb;
        tT(w_count + i) = T(first);
        tb(w_count + i) = b(first);
        tbetta(w_count + i) = betta(first);
        first = next + 1;
    end
end

end