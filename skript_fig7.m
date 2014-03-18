% !!! Calculation of parameter maps takes some time !!! 
%% calculate parameter maps for all 6 volunteer data sets (approx 35 min).

%for each volunteer
for p = 1:6
load(['data/volunteer' num2str(p)]);

%3x3 median filter
for i = 1:size(data, 1)
    data(i, :, :) = medfilt2(reshape(data(i, :), size(data, 2), size(data, 3)), [3 3]);
end

%average over diffusion directions
[data, betta, b, T] = trace(data, betta, b, T);
b0 = mean(data(b == 0, mask > 0), 1);
data = data(b > 0, mask > 0);
T = T(b > 0);
betta = betta(b > 0);
b = b(b > 0);

%normalize to average b0 image
for i = 1:size(data, 1)
    data(i, :) = data(i, :) ./ b0;  
end

% load phase distributions and set acquisition parameters
global generic
generic = load('phd/generic_v2.mat');
generic.('T') = T;
generic.('b') = b;
generic.('betta') = betta;
generic.('Db') = ones(size(b)) * 1.6;

%initialize variables for pixelwise fit results
Dhighb = zeros(1, size(data, 2));
fhighb = zeros(1, size(data, 2));
D = zeros(1, size(data, 2));
f = zeros(1, size(data, 2));
Dstar = zeros(1, size(data, 2));
tau = zeros(1, size(data, 2));
r = zeros(1, size(data, 2));
v = zeros(1, size(data, 2));

% selector for high b values and bipolar gradients
biphighb = ((b >= 200) & (betta == 1));
options = optimset('Display', 'off');

tic
for i = 1:size(data, 2)
    
    %[param(1:2), ~] = fminsearch(@(x) norm(data(biphighb, i) - (1-x(2))*exp(-b(biphighb) * x(1) / 1000))^2, [0.2, 2], options);
    lb = [0, 0];
    ub = [4, 1];
    param(1:2) = lsqnonlin(@(x) (data(biphighb, i) - (1-x(2))*exp(-b(biphighb) * x(1) / 1000)), [0.2, 2], lb, ub, options);
    Dhighb(i) = param(1);
    fhighb(i) = param(2);
    %[param(3:6), param(7)] = fminsearch(@(x) norm(data(:, i) - generic_model(x))^2, [param(1:2), 2, 4], options);
    lb = [0, 0, 1/100, 0];
    ub = [4, 1, 1599/100, 1000];
    [param(3:6), param(7)] = lsqnonlin(@(x) data(:, i) - generic_model(x), [param(1:2), 2, 4], lb, ub, options);
    D(i) = param(3);
    f(i) = param(4);
    tau(i) = param(5) * 100;
    v(i) = param(6);
    Dstar(i) = tau(i) *  v(i) * v(i) / 6;
    r(i) = param(7);
    
    waitbar(i / size(data, 2))
end
toc

temp = zeros(size(mask));
temp(mask > 0) = b0;
b0 = temp;
temp(mask > 0) = D;
D = temp;
temp(mask > 0) = f;
f = temp;
temp(mask > 0) = r;
r = temp;
temp(mask > 0) = v;
v = temp;
temp(mask > 0) = Dstar;
Dstar = temp;
temp(mask > 0) = tau;
tau = temp;
temp(mask > 0) = Dhighb;
Dhighb = temp;
temp(mask > 0) = fhighb;
fhighb = temp;
clear lb_highb ub_highb lb ub options1 highbdata highb

if (~exist('maps', 'dir'))
    mkdir('maps');
end
save(['maps/maps_volunteer' num2str(p)], 'b0', 'D', 'f', 'tau', 'v', 'Dstar', 'Dhighb', 'fhighb', 'r');

end

%% plot color maps for all probands into figure:
figure(1)
ha = tight_subplot(6, 6, 0, 0, 0);
Dmax = 4;
fmax = 1;
taumax = 600;
Dstarmax = 2000;
vmax = 12.5;
rmax = 5;
cmap = jet(128);
cmap2 = gray(128);
cmap = [flipud(cmap2); cmap];
colormap(cmap);
set(gcf, 'Color', 'k')

selecty1 = [27; 20; 22; 22; 18; 22];
selecty2 = [75; 65; 64; 63; 61; 61];
selectx1 = [15; 12; 18; 11; 15; 15];
selectx2 = [88; 93; 85; 93; 89; 85];

for i = 1:6

selecty = selecty1(i):selecty2(i);
selectx = selectx1(i):selectx2(i);
maps = load(['maps/maps_volunteer' num2str(i) '.mat']);
maps.('D')(maps.('b0') == 0) = -0.01;
maps.('f')(maps.('b0') == 0) = -0.01;
maps.('tau')(maps.('b0') == 0) = -0.01;
maps.('Dstar')(maps.('b0') == 0') = -0.01;
maps.('v')(maps.('b0') == 0) = -0.01;
maps.('b0')(maps.('b0') == 0) = 0.01;

axes(ha(i));
imshow(-maps.('b0')(selecty, selectx), [-1 1] * max(max(maps.('b0')(selecty, selectx))), 'Colormap', cmap);
axis off
axes(ha(i + 6));
imshow(maps.('D')(selecty, selectx), [-1 1] * Dmax, 'Colormap', cmap);
axis off
axes(ha(i + 12));
imshow(maps.('f')(selecty, selectx), [-1 1] * fmax, 'Colormap', cmap);
axis off
axes(ha(i + 18));
imshow(maps.('tau')(selecty, selectx), [-1 1] * taumax, 'Colormap', cmap);
axis off
axes(ha(i + 24));
imshow(maps.('v')(selecty, selectx), [-1 1] * vmax, 'Colormap', cmap);
axis off
axes(ha(i + 30));
imshow(maps.('Dstar')(selecty, selectx), [-1 1] * Dstarmax, 'Colormap', cmap);
axis off

end
