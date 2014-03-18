
%% display signal attenuation in ROIs in pancreas and liver and fit the
% generic IVIM model to the data

% initialize variables for averaged ROI values
betta_avg = reshape(ones(20, 1) * [2 1], 1, 40)';
T_avg = reshape([ones(4, 1) * [40 70 70 100 100] ones(4, 1) * [40 70 70 100 100]], 1, 40)';
bs = [10 25 50 100 200 300 400 500];
b_avg = [bs(1:4) bs bs bs(1:4) bs bs]';
roi_data = zeros(numel(b_avg), 6);
pancreas =  zeros(numel(b_avg), 1);
s_pancreas =  zeros(numel(b_avg), 1); %for error bars
liver =  zeros(numel(b_avg), 1);
s_liver =  zeros(numel(b_avg), 1);    %for error bars

% load pancreas ROI data and prepare averaging over volunteers:
roi = 1; %1=pancreas, 2=liver
for v = 1:6
    load(['data/volunteer' num2str(v) '.mat']);
    [data, betta, b, T] = trace(data, betta, b, T);
    data = mean(data(:, rois == roi), 2);
    data = data(b > 0) ./ data(floor((find(b > 0) - 1) / 3) * 3 + 1);
    betta = betta(b > 0);
    T = T(b > 0);
    b = b(b > 0);
    % in case of a particular acquisition parameter set was acquired
    % several times for a volunteer, average: 
    for i = 1:numel(b_avg)
        if (numel(find((b == b_avg(i) & T == T_avg(i) & betta == betta_avg(i)))))
            roi_data(i, v) = mean(data(b == b_avg(i) & T == T_avg(i) & betta == betta_avg(i)));
        end
    end
end

% average over volunteers while taking into account that a particular 
% set of acquisition parameters might have been forgotten during
% acquisition (default value 0 not overwritten)
for i = 1:numel(b_avg)
    pancreas(i) = mean(roi_data(i, roi_data(i, :) > 0));
    s_pancreas(i) = std(roi_data(i, roi_data(i, :) > 0), 0, 2);
end

% same for liver ....

% load liver ROI data and average over volunteers:
roi = 2; %1=pancreas, 2=liver
for v = 1:6
    load(['data/volunteer' num2str(v) '.mat']);
    [data, betta, b, T] = trace(data, betta, b, T);
    data = mean(data(:, rois == roi), 2);
    data = data(b > 0) ./ data(floor((find(b > 0) - 1) / 3) * 3 + 1);
    betta = betta(b > 0);
    T = T(b > 0);
    b = b(b > 0);
    for i = 1:numel(b_avg)
        if (numel(find((b == b_avg(i) & T == T_avg(i) & betta == betta_avg(i)))))
            roi_data(i, v) = mean(data(b == b_avg(i) & T == T_avg(i) & betta == betta_avg(i)));
        end
    end
end

for i = 1:numel(b_avg)
    liver(i) = mean(roi_data(i, roi_data(i, :) > 0));
    s_liver(i) = std(roi_data(i, roi_data(i, :) > 0), 0, 2);
end

b = b_avg;
T = T_avg;
betta = betta_avg;
clear data b_avg T_avg betta_avg roi_data i v rois roi mask bs

%plot signal attenuation in ROIs with error bars:
FC40 = betta == 2 & T == 40;
FC70 = betta == 2 & T == 70;
FC100 = betta == 2 & T == 100;
BP40 = betta == 1 & T == 40;
BP70 = betta == 1 & T == 70;
BP100 = betta == 1 & T == 100;

figure(1)
hold off
errorbar(b(FC40), liver(FC40), s_liver(FC40), 'ro');
hold on
errorbar(b(FC70), liver(FC70), s_liver(FC70), 'g^');
errorbar(b(FC100), liver(FC100), s_liver(FC100), 'ks');
errorbar(b(BP40), liver(BP40), s_liver(BP40), 'rh');
errorbar(b(BP70), liver(BP70), s_liver(BP70), 'gp');
errorbar(b(BP100), liver(BP100), s_liver(BP100), 'kd');
xlabel('b (s/mm²)');
ylabel('S(b) / S(0)');
title('Signal attenuation for liver');
legend('FC40', 'FC70', 'FC100', 'BP40', 'BP70', 'BP100');

figure(2)
hold off
errorbar(b(FC40), pancreas(FC40), s_pancreas(FC40), 'ro');
hold on
errorbar(b(FC70), pancreas(FC70), s_pancreas(FC70), 'g^');
errorbar(b(FC100), pancreas(FC100), s_pancreas(FC100), 'ks');
errorbar(b(BP40), pancreas(BP40), s_pancreas(BP40), 'rh');
errorbar(b(BP70), pancreas(BP70), s_pancreas(BP70), 'gp');
errorbar(b(BP100), pancreas(BP100), s_pancreas(BP100), 'kd');
xlabel('b (s/mm²)');
ylabel('S(b) / S(0)');
title('Signal attenuation for pancreas');
legend('FC40', 'FC70', 'FC100', 'BP40', 'BP70', 'BP100');


%perform fitting for liver data:
global generic
generic = load('phd/generic_v2.mat');
generic.('T') = T;
generic.('b') = b;
generic.('betta') = betta;
generic.('Db') = ones(size(b)) * 1.6;

biphighb = ((b >= 200) & (betta == 1));
options = optimset('Display', 'off');
[param(1:2), param(3)] = fminsearch(@(x) norm(liver(biphighb) - (1-x(2))*exp(-b(biphighb) * x(1) / 1000))^2, [0.2, 2], options);
[param(4:7), param(8)] = fminsearch(@(x) norm(liver - generic_model(x))^2, [param(1:2), 2, 4]);

resi = param(8) / (numel(generic.('b')) - 4);
jacky = jacobianest(@(x) (liver - generic_model(x)), param(4:7));
param(9:12) = sqrt(diag(resi * inv(jacky'*jacky))); %#ok<MINV>

fprintf('Best fit liver: D=%f±%f, f=%f±%f, tau=%f±%f, v=%f±%f\n', param(4), ...
    param(9), param(5), param(10), param(6) * 100, param(11) * 100, param(7), param(12));
%Best fit liver: D=1.215519±0.040583, f=0.346731±0.007947, tau=143.783913±9.717227, v=4.595080±0.336984

%alternative code for lsqnonlin / nlparci
%lb = [0, 0];
%ub = [4, 1];
%[param(1:2), param(3), ~, ~, ~] = lsqnonlin(@(x) (liver(biphighb) - (1-x(2))*exp(-b(biphighb) * x(1) / 1000)), [0.2, 2], lb, ub, options);
%lb = [0, 0, 1/100, 0];
%ub = [4, 1, 1599/100, 1000];
%[param(4:7), param(8), resi, ~, ~, ~, jacky] = lsqnonlin(@(x) liver - generic_model(x), [param(1:2), 2, 4], lb, ub, options);
%ci = nlparci(param(4:7), resi, jacky, 0.05);
%fprintf('Best fit liver: D=%f±%f, f=%f±%f, tau=%f±%f, v=%f±%f\n', param(4), ...
%    (ci(1,2) - ci(1,1)) / 2, param(5), (ci(2,2) - ci(2,1)) / 2, param(6) * 100, ...
%    (ci(3,2) - ci(3,1)) / 2 * 100, param(7), (ci(4,2) - ci(4,1)) / 2);

%plot fit solution, acquisition parameters need to be set in fields of
%global variable generic:
generic.('b') = (0:550)';
generic.('Db') = ones(size(generic.('b'))) * 1.6;
generic.('betta') = ones(size(generic.('b'))) * 2;

generic.('T') = ones(size(generic.('b'))) * 40;
FC40_fit = generic_model(param(4:7));
generic.('T') = ones(size(generic.('b'))) * 70;
FC70_fit = generic_model(param(4:7));
generic.('T') = ones(size(generic.('b'))) * 100;
FC100_fit = generic_model(param(4:7));
generic.('betta') = ones(size(generic.('b')));
MP100_fit = generic_model(param(4:7));
generic.('T') = ones(size(generic.('b'))) * 70;
MP70_fit = generic_model(param(4:7));
generic.('T') = ones(size(generic.('b'))) * 40;
MP40_fit = generic_model(param(4:7));

figure(1)
plot(FC40_fit, 'r');
plot(FC70_fit, 'g');
plot(FC100_fit, 'k');
plot(MP40_fit, 'r');
plot(MP70_fit, 'g');
plot(MP100_fit, 'k');

%perform fitting for pancreas data:
generic = load('phd/generic_v2.mat');
generic.('T') = T;
generic.('b') = b;
generic.('betta') = betta;
generic.('Db') = ones(size(b)) * 1.6;

biphighb = ((b >= 200) & (betta == 1));
options = optimset('Display', 'off');
[param(1:2), param(3)] = fminsearch(@(x) norm(pancreas(biphighb) - (1-x(2))*exp(-b(biphighb) * x(1) / 1000))^2, [0.2, 2], options);
[param(4:7), param(8)] = fminsearch(@(x) norm(pancreas - generic_model(x))^2, [param(1:2), 2, 4]);

resi = param(8) / (numel(generic.('b')) - 4);
jacky = jacobianest(@(x) (pancreas - generic_model(x)), param(4:7));
param(9:12) = sqrt(diag(resi * inv(jacky'*jacky))); %#ok<MINV>

fprintf('Best fit pancreas: D=%f±%f, f=%f±%f, tau=%f±%f, v=%f±%f\n', param(4), ...
    param(9), param(5), param(10), param(6) * 100, param(11) * 100, param(7), param(12));
%Best fit pancreas: D=1.876034±0.080872, f=0.287270±0.015398, tau=224.136403±47.215061, v=3.911401±0.544393

%alternative code for lsqnonlin / nlparci
%lb = [0, 0];
%ub = [4, 1];
%[param(1:2), param(3), ~, ~, ~] = lsqnonlin(@(x) (pancreas(biphighb) - (1-x(2))*exp(-b(biphighb) * x(1) / 1000)), [0.2, 2], lb, ub, options);
%lb = [0, 0, 1/100, 0];
%ub = [4, 1, 1599/100, 1000];
%[param(4:7), param(8), resi, ~, ~, ~, jacky] = lsqnonlin(@(x) pancreas - generic_model(x), [param(1:2), 2, 4], lb, ub, options);
%ci = nlparci(param(4:7), resi, jacky, 0.05);
%fprintf('Best fit pancreas: D=%f±%f, f=%f±%f, tau=%f±%f, v=%f±%f\n', param(4), ...
%    (ci(1,2) - ci(1,1)) / 2, param(5), (ci(2,2) - ci(2,1)) / 2, param(6) * 100, ...
%    (ci(3,2) - ci(3,1)) / 2 * 100, param(7), (ci(4,2) - ci(4,1)) / 2);

generic.('b') = (0:550)';
generic.('Db') = ones(size(generic.('b'))) * 1.6;
generic.('betta') = ones(size(generic.('b'))) * 2;

generic.('T') = ones(size(generic.('b'))) * 40;
FC40_fit = generic_model(param(4:7));
generic.('T') = ones(size(generic.('b'))) * 70;
FC70_fit = generic_model(param(4:7));
generic.('T') = ones(size(generic.('b'))) * 100;
FC100_fit = generic_model(param(4:7));
generic.('betta') = ones(size(generic.('b')));
MP100_fit = generic_model(param(4:7));
generic.('T') = ones(size(generic.('b'))) * 70;
MP70_fit = generic_model(param(4:7));
generic.('T') = ones(size(generic.('b'))) * 40;
MP40_fit = generic_model(param(4:7));

figure(2)
plot(FC40_fit, 'r');
plot(FC70_fit, 'g');
plot(FC100_fit, 'k');
plot(MP40_fit, 'r');
plot(MP70_fit, 'g');
plot(MP100_fit, 'k');

