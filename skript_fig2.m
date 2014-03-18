
%comparison of GPA and exact signal attenuation in case of IVIM for flow 
%compensated and bipolar gradients for a single velocity and a parabolic 
%flow profile

%simulation parameters:
T = 70;        %ms
tau = 180;      %ms
v = 4.25;         %mm/s
b = 0:550;     %s/mm²
options = optimset('Display', 'off');

%print nominal value of Dstar (µm²/ms):
fprintf('Dstar_nominal = %f\n', tau * v * v / 6);

%load phase distributions for a single velocity component:
temp = load('phd/generic.mat', 'phis', 'fc1', 'bip', 'N');
N = temp.('N');
phis = temp.('phis');

%obtain phase distribution for specified T/tau ratio via linear interpolation:
Ns = T / tau;
iceil = max(0, ceil(log(Ns/N(2)) / log(max(N)/N(2)) * (numel(N)-2))) + 2;
ifloor = max(-1, floor(log(Ns/N(2)) / log(max(N)/N(2)) * (numel(N)-2))) + 2;
c = (Ns - N(ifloor)) ./ (N(iceil) - N(ifloor)) * ones(size(phis));
c(N(iceil) == N(ifloor), :) = 1;
bip = (1-c) .* temp.('bip')(ifloor, :) + c .* temp.('bip')(iceil, :);
fc1 = (1-c) .* temp.('fc1')(ifloor, :) + c .* temp.('fc1')(iceil, :);

%calculate signal attenuation with and without Gaussian phase approximation

%GPA:
F_fc1 = exp(-b * T * v * v / 1000 * sum(phis .^2 .* fc1) / 2);
F_bip = exp(-b * T * v * v / 1000 * sum(phis .^2 .* bip) / 2);

%exact:
F_fc1_cor = abs(sum((ones(size(b')) * fc1) .* cos(sqrt(b' * T * v * v / 1000) * phis), 2));
F_bip_cor = abs(sum((ones(size(b')) * bip) .* cos(sqrt(b' * T * v * v / 1000) * phis), 2));

%print apparent Dstar for GPA signal attenuations
fprintf('Dstar_fc1 = %f\n', lsqnonlin(@(x) exp(-b * x / 1000) - F_fc1, tau * v * v / 6, 0, 10000, options));
fprintf('Dstar_bip = %f\n', lsqnonlin(@(x) exp(-b * x / 1000) - F_bip, tau * v * v / 6, 0, 10000, options));

%plot attenuation curves for single velocity:
figure(1)
hold off
plot(b, F_fc1, ':r');
hold on
plot(b, F_fc1_cor, 'r');
plot(b, F_bip, ':b');
plot(b, F_bip_cor, 'b');
title('IVIM signal for a single velocity with and without GPA (\tau=180ms, T=70ms, v=4.25mm/s)')
legend('F_{fc1} (GPA)', 'F_{fc1} (exact)', 'F_{bip} (GPA)', 'F_{bip} (exact)')
xlabel('b (s/mm²)');
ylabel('F(b, T, \tau, v)');

export_v1 = [F_fc1' F_fc1_cor F_bip' F_bip_cor b'];

%analogous procedure for parabolic flow profile:

%load phase distributions with integrated parabolic velocity profile
temp = load('phd/generic_v2.mat', 'phis', 'fc1', 'bip', 'N');
N = temp.('N');
phis = temp.('phis');

%obtain phase distribution for specified T/tau ratio:
Ns = T / tau;
iceil = max(0, ceil(log(Ns/N(2)) / log(max(N)/N(2)) * (numel(N)-2))) + 2;
ifloor = max(-1, floor(log(Ns/N(2)) / log(max(N)/N(2)) * (numel(N)-2))) + 2;
c = (Ns - N(ifloor)) ./ (N(iceil) - N(ifloor)) * ones(size(phis));
c(N(iceil) == N(ifloor), :) = 1;
bip = (1-c) .* temp.('bip')(ifloor, :) + c .* temp.('bip')(iceil, :);
fc1 = (1-c) .* temp.('fc1')(ifloor, :) + c .* temp.('fc1')(iceil, :);

%calculate signal attenuation with and without Gaussian phase approximation
F_fc1 = exp(-b * T * v * v / 1000 * sum(phis .^2 .* fc1) / 2);
F_bip = exp(-b * T * v * v / 1000 * sum(phis .^2 .* bip) / 2);
F_fc1_cor = abs(sum((ones(size(b')) * fc1) .* cos(sqrt(b' * T * v * v / 1000) * phis), 2));
F_bip_cor = abs(sum((ones(size(b')) * bip) .* cos(sqrt(b' * T * v * v / 1000) * phis), 2));

%plot results
figure(2)
hold off
plot(b, F_fc1, ':r');
hold on
plot(b, F_fc1_cor, 'r');
plot(b, F_bip, ':b');
plot(b, F_bip_cor, 'b');
title('IVIM signal for a parabolic flow profile with and without GPA (\tau=180ms, T=70ms, v=4.25mm/s)')
legend('F_{fc1} (GPA)', 'F_{fc1} (exact)', 'F_{bip} (GPA)', 'F_{bip} (exact)')
xlabel('b (s/mm²)');
ylabel('F(b, T, \tau, v)');

export_v2 = [F_fc1' F_fc1_cor F_bip' F_bip_cor b'];

%display apparent Dstar
fprintf('Dstar_fc1_v2 = %f\n', lsqnonlin(@(x) exp(-b * x / 1000) - F_fc1, tau * v * v / 6, 0, 10000, options));
fprintf('Dstar_bip_v2 = %f\n', lsqnonlin(@(x) exp(-b * x / 1000) - F_bip, tau * v * v / 6, 0, 10000, options));