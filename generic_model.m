% calculates the signal attenuation in the IVIM model for the model
% parameters specified in x and the acquisition parameters set via the
% global variable "generic", which also has to contain the phase 
% distributions for the bipolar ('bip') and flow compensated ('fc1') 
% gradient profiles.
%
% model parameters are specified as follows:
%   D   = x(1)       µm²/ms
%   f   = x(2)
%   tau = x(3) * 100 -> must be specified in x in multiples of 0.1 s
%   v   = x(4)       mm/s
%
% The output is of the same length as the fields specifying the acquisition
% parameters. This version is optimized for speeding up the signal 
% calculation for the used imaging parameters such that only diffusion 
% times T=40ms, T=70ms and T=100ms are considered. By implementing a loop 
% over generic.('T') arbitrary values for T can be evaluated at the expense
% of computational time. The fields that need to be set in "generic" are:
%
% The acquisition parameters need to be of the same size and determine the
% output size:
% - generic.('T'): diffusion times
% - generic.('b'): b-values
% - generic.('betta'): bipolar (bip) or flow compensated (fc1)
%
% The phase distributions are stored in the following set of fields:
% - generic.('N'): ratios N=T/tau for which phase distributions are stored
% - generic.('phis'): centers of histogram bins of the phase distributions
% - generic.('bip'): bipolar 
% - generic.('fc1'): flow compensated 
%
function S = generic_model(x)

global generic

D = x(1);
f = x(2);
tau = abs(x(3)) * 100; % ms
v = x(4); % mm/s

phds = zeros(length(generic.('T')), length(generic.('phis')));
None = ones(size(generic.('T')));

xmax = max(generic.('phis'));

% First handle acquisition parameters with T = 100 ms.
iceil = find(generic.('N') >= 100/tau, 1);
% if N=T/tau is larger than all values of N for which phase distributions
% are available (default N_max = 100), then use Gaussian approximation:
if (numel(iceil) == 0)
    logos = (generic.('T') == 100);
    phds(logos, :) = None(logos) * ...
        (erf((min(xmax, generic.('phis') + generic.('phis')(2) / 2)) * sqrt(150 / tau)) ...
        - erf(max(0, generic.('phis') - generic.('phis')(2) / 2) * sqrt(150 / tau)));
    
% otherwise do a linear interpolation between the two phase distributions
% corresponding to the two surrounding values of N.
else
    ifloor = find(generic.('N') < 100/tau, 1, 'last');
    c = (100/tau - generic.('N')(ifloor)) / (generic.('N')(iceil) - generic.('N')(ifloor));
    logos = (generic.('T') == 100) & (generic.('betta') == 1);
    phds(logos, :) = None(logos) * ...
        ((1-c) * generic.('bip')(ifloor, :) + c * generic.('bip')(iceil, :));
    logos = (generic.('T') == 100) & (generic.('betta') == 2);
    phds(logos, :) = None(logos) * ...
        ((1-c) * generic.('fc1')(ifloor, :) + c * generic.('fc1')(iceil, :));
end

% Then handle acquisition parameters with T = 70 ms (basically same as
% above)...
iceil = find(generic.('N') >= 70/tau, 1);
if (numel(iceil) == 0)
    logos = (generic.('T') == 70);
    phds(logos, :) = None(logos) * ...
        (erf((min(xmax, generic.('phis') + generic.('phis')(2) / 2)) * sqrt(105 / tau)) ...
        - erf(max(0, generic.('phis') - generic.('phis')(2) / 2) * sqrt(105 / tau)));
else
    ifloor = find(generic.('N') < 70/tau, 1, 'last');
    c = (70/tau - generic.('N')(ifloor)) / (generic.('N')(iceil) - generic.('N')(ifloor));
    logos = (generic.('T') == 70) & (generic.('betta') == 1);
    phds(logos, :) = None(logos) * ...
        ((1-c) * generic.('bip')(ifloor, :) + c * generic.('bip')(iceil, :));
    logos = (generic.('T') == 70) & (generic.('betta') == 2);
    phds(logos, :) = None(logos) * ...
        ((1-c) * generic.('fc1')(ifloor, :) + c * generic.('fc1')(iceil, :));
end

% Finally T = 40 ms
iceil = find(generic.('N') >= 40/tau, 1);
if (numel(iceil) == 0)
    logos = (generic.('T') == 40);
    phds(logos, :) = None(logos) * ...
        (erf((min(xmax, generic.('phis') + generic.('phis')(2) / 2)) * sqrt(60 / tau)) ...
        - erf(max(0, generic.('phis') - generic.('phis')(2) / 2) * sqrt(60 / tau)));
else
    ifloor = find(generic.('N') < 40/tau, 1, 'last');
    c = (40/tau - generic.('N')(ifloor)) / (generic.('N')(iceil) - generic.('N')(ifloor));
    logos = (generic.('T') == 40) & (generic.('betta') == 1);
    phds(logos, :) = None(logos) * ...
        ((1-c) * generic.('bip')(ifloor, :) + c * generic.('bip')(iceil, :));
    logos = (generic.('T') == 40) & (generic.('betta') == 2);
    phds(logos, :) = None(logos) * ...
        ((1-c) * generic.('fc1')(ifloor, :) + c * generic.('fc1')(iceil, :));
end

% Calculate signal attenuation of the perfusion fraction from the phase distributions:
F = abs(sum(phds .* cos(v * sqrt(generic.('b').*generic.('T') / 1000) * generic.('phis')), 2));

% Calculate signal in the two-compartment model:
S = ((1 - f) * exp(-generic.('b') * D / 1000) + f * exp(-generic.('b') .* generic.('Db') / 1000) .* F);
        
end
