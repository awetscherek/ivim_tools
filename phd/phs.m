function retval = phs(x)
% matrix version of the "Pseudo"-HeaviSide function

    retval = floor(1 ./ (1 + exp(-10000000 * x)));

end
