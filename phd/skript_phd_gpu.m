%          !!! This script must be executed from within !!!  
%          !!! the folder containing the phd_*.m files. !!!

%contents:
%    1. skript for calculating phase histograms 
%    2. normalization of phase histograms
%    3. use exact solution for N = 1
%    4. apply parabolic velocity profile

%% 1. calculation of phase histograms for different N (GPU version)
% total run time approx. 7h 45 min on i7-3820 / GTX 660 Ti

% settings:
 % select weighting function for gradient profile 
profiles = ['bip';'fc0';'fc1']; % bip = bipolar profile
                                % fc0 = simplest flow comp waveform (FEER)
                                % fc1 = flow comp waveform used in 
                                %       manuscript.
seed     =     220784; % reproducibility.
minN     =          1; % smallest N (typically 1).
maxN     =        100; % largest N (>~30).
Nsteps   =       1024; % number of different N.
phisteps =       1024; % number of phase histogram columns.
paths    =   64000000; % number of particle paths.
output   = 'base.mat'; % file to store phase distributions.

%calculate different Ns
Ns = exp((0:(Nsteps-1)) / (Nsteps-1) * log(maxN / minN)) * minN; 

%column edges for phase histograms
phis = sqrt(3) * ((0:phisteps) - 0.5) ... %-0.5 because of histc
       / (phisteps-1) / 2;

%initialize random number generator on GPU:
sg = parallel.gpu.RandStream('CombRecursive', 'Seed', seed);
parallel.gpu.RandStream.setGlobalStream(sg);

clear minN maxN Nsteps phisteps

% for each profile:
for p = 1:size(profiles, 1)
    
    %total calculation time counter for current profile:
    T = 0;
    
    %initialize variable for phase histograms
    eval([profiles(p, :) ' = zeros(numel(Ns), numel(phis));']);

    %choose path weightings
    eval(['cfun = @c_' profiles(p, :) ';'] );

    for Nind = numel(Ns):-1:1

        tic;
        %temporary variable for particle phases and random starting points:
        phases = gpuArray.zeros(1, paths, 'single');
        r = gpuArray.rand(1, paths, 'single');
        %for each segment:
        for k = 0:ceil(Ns(Nind))
            c = arrayfun(cfun, r + k, Ns(Nind));
            phases = phases + arrayfun(@(x,y) x.*(1-2*y), ...
                     c, gpuArray.rand(1, paths, 'single'));
            clear c
        end
        eval([profiles(p, :) '(Nind, :) = ' profiles(p, :) ...
            '(Nind, :) + histc(gather(abs(phases)), phis);']);
        clear r phases
        T = T + toc;
    end
    
    clear Nind T cfun k
    
    eval([profiles(p, :) ' = ' profiles(p, :) '(:, 1:(numel(phis)-1));']);
    
    fprintf('total time for profile %s: %f s\n', profiles(p, :), T);
    % total time for profile bip:  7689.670926 s
    % total time for profile fc0:  7935.187316 s
    % total time for profile fc1: 12353.691162 s

end

% save centers of histogram columns
phis = sqrt(3) * (0:(numel(phis)-2)) / (numel(phis)-2) / 2;
N = Ns';
save(output, 'N', 'phis', 'profiles');
for p = 1:size(profiles, 1)
    save(output, profiles(p, :), '-append');
end
clear p Ns output phis N profiles sg

%% 2. normalize phase distributions 
input = 'base.mat';
output = 'generic.mat';

phd_normalize(input, output);
clear input output

%% 3. use exact solution for N = 1
input = 'generic.mat';
output = 'generic_N1.mat';

phd_exact_N1(input, output);
clear input output

%% 4. apply parabolic velocity profile = constant velocity distribution
% ~45 min on i7-3820 / GTX 660 Ti
input = 'generic_N1.mat';
output = 'generic_v2.mat';

% phd for N in range between 0 and 1 need to be saved independently, when
% a velocity distribution is assumed.
minN = 0.01;
maxN = 100;
Nsteps = 2047;

vsteps = 2048; % discretization steps

%calculate different Ns
N2 = [0 exp((0:(Nsteps-1)) / (Nsteps-1) * log(maxN / minN)) * minN]'; 
tic
phd_apply_v2(input, output, N2, vsteps);
toc
clear input output minN maxN Nsteps vsteps N2
