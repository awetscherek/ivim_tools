
%% load roi data for pipe phantom with and without flow

%only evaluate pipe-ROI
roi = 2; %1=static bottle, 2=flow tube

%load data without flow
load('data/phantom_static.mat');
[data_static, betta_static, b_static, T_static] = trace14(data, betta, b, T);
tube_static = mean(data_static(:, rois == roi), 2);
tube_static = tube_static(b_static > 0) / mean(tube_static(b_static == 0));
betta_static = betta_static(b_static > 0);
T_static = T_static(b_static > 0);
b_static = b_static(b_static > 0);

%load data with flow
load('data/phantom_flow.mat');
[data_flow, betta_flow, b_flow, T_flow] = trace14(data, betta, b, T);
tube_flow = mean(data_flow(:, rois == roi), 2);
tube_flow = tube_flow(b_flow > 0) / mean(tube_flow(b_flow == 0));
betta_flow = betta_flow(b_flow > 0);
T_flow = T_flow(b_flow > 0);
b_flow = b_flow(b_flow > 0);
clear roi rois data betta b T mask data_static data_flow

%plot flow comp and bipolar data for flow into figure:
figure(1)
hold off
plot(b_flow((betta_flow == 2) & (T_flow == 40)), tube_flow((betta_flow == 2) & (T_flow == 40)), 'ro');
hold on
plot(b_flow((betta_flow == 2) & (T_flow == 70)), tube_flow((betta_flow == 2) & (T_flow == 70)), 'g^');
plot(b_flow((betta_flow == 2) & (T_flow == 100)), tube_flow((betta_flow == 2) & (T_flow == 100)), 'ks');
plot(b_flow((betta_flow == 1) & (T_flow == 40)), tube_flow((betta_flow == 1) & (T_flow == 40)), 'rs');
plot(b_flow((betta_flow == 1) & (T_flow == 70)), tube_flow((betta_flow == 1) & (T_flow == 70)), 'g<');
plot(b_flow((betta_flow == 1) & (T_flow == 100)), tube_flow((betta_flow == 1) & (T_flow == 100)), 'k>');

%plot bipolar data without flow into figure:
plot(b_static((betta_static == 1) & (T_static == 40)), tube_static((betta_static == 1) & (T_static == 40)), 'mo');
plot(b_static((betta_static == 1) & (T_static == 70)), tube_static((betta_static == 1) & (T_static == 70)), 'b<');
plot(b_static((betta_static == 1) & (T_static == 100)), tube_static((betta_static == 1) & (T_static == 100)), 'c>');

legend('FC40 flowing', 'FC70 flowing', 'FC100 flowing', 'BP40 flowing', ...
    'BP70 flowing', 'BP100 flowing', 'BP40 static', 'BP70 static', 'BP100 static');
title('Signal attenuation measured in pipe phantom with and without flow');
xlabel('b (s/mm²)');
ylabel('S(b)');

%prepare fit of a monoexponential decay curve:
b = [b_flow(betta_flow == 2)' b_static(betta_static == 1)'];
S = [tube_flow(betta_flow == 2)' tube_static(betta_static == 1)'];

options = optimset('Display', 'off');
D_fit = fminsearch(@(x) norm(exp(-b * x / 1000) - S)^2, 2, options);

%print fitted diffusion coefficient and plot associated signal attenuation
%curve:
fprintf('D_fit = %f\n', D_fit);
plot(0:500, exp(-(0:500)/1000 * D_fit), 'k:');
