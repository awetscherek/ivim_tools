
%% displays phase distributions for selected values of Ns
% adds analytic solutions for N=1 and the pseudo-diffusion limit for the
% largest value in Ns

% selected values
Ns = [1 3 10 30];

% load phase distribtions
temp = load('phd/generic.mat', 'fc1', 'bip', 'N', 'phis');
N = temp.('N');
phis = temp.('phis');

% calculate closest values of N for which phase distributions are known:
iceil = max(0, ceil(log(Ns/N(2)) / log(max(N)/N(2)) * (numel(N)-2))) + 2;
ifloor = max(-1, floor(log(Ns/N(2)) / log(max(N)/N(2)) * (numel(N)-2))) + 2;

% determine weighting:
c = (Ns' - N(ifloor)) ./ (N(iceil) - N(ifloor)) * ones(size(phis));
c(N(iceil) == N(ifloor), :) = 1;

% calculate phase distributions as a linear combination:
bip = (1-c) .* temp.('bip')(ifloor, :) + c .* temp.('bip')(iceil, :);
fc1 = (1-c) .* temp.('fc1')(ifloor, :) + c .* temp.('fc1')(iceil, :);

%account for histogram bin for phi=0 is only of half width:
bip(:, 1) = bip(:, 1) * 2;
fc1(:, 1) = fc1(:, 1) * 2;

addpath('phd');

figure(1)
hold off
%plot flow compensated phase distributions in figure 1
plot(phis, circshift(fc1', [0 -1]) / sqrt(3) * (numel(phis)-1), 'LineWidth', 2)
hold on
%plot analytic solutions:
plot(phis, phd_fc1_N1(phis), ':k', 'LineWidth', 2)
plot(phis, phd_lim_Nhigh(phis, max(Ns)), '-.k', 'LineWidth', 2)
axis([0 sqrt(3)/2 0 4]);
set(gca, 'XDir', 'reverse');
legend('N=3 sim', 'N=10 sim', 'N=30 sim', 'N=1 sim', 'N=1 exact', 'N=30 limit', 'Location', 'NorthWest');
xlabel('\phi');
ylabel('\rho(\phi, N)');
title('Flow compensated phase distributions for selected N');

%same for bipolar phase distributions
figure(2)
hold off
plot(phis, circshift(bip', [0 -1]) / sqrt(3) * (numel(phis)-1), 'LineWidth', 2)
hold on
plot(phis, phd_bip_N1(phis), ':k', 'LineWidth', 2)
plot(phis, phd_lim_Nhigh(phis, max(Ns)), '-.k', 'LineWidth', 2)
axis([0 sqrt(3)/2 0 4]);
xlabel('\phi');
ylabel('\rho(\phi, N)');
title('Bipolar phase distributions for selected N');
legend('N=3 sim', 'N=10 sim', 'N=30 sim', 'N=1 sim', 'N=1 exact', 'N=30 limit', 'Location', 'NorthEast');

bip = bip ./ sqrt(3) * (numel(phis)-1);
fc1 = fc1 ./ sqrt(3) * (numel(phis)-1);

export = [bip(1, :)' fc1(1, :)' bip(2, :)' fc1(2, :)' bip(3, :)' fc1(3, :)' ...
          bip(4, :)' fc1(4, :)'  phd_lim_Nhigh(phis', max(Ns)) ...
          phd_bip_N1(phis')  phd_fc1_N1(phis')];   
export = [-flipud(phis(2:numel(phis))') flipud(export(2:size(export,1), :)) ; phis' export];

rmpath('phd');