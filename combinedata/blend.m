%% Read in data
close all
clear
addpath(".\Data\");
%prettygraphs;
filename = "kelley_627_multi2";
rawdata  = csvread(filename + ".csv");
rawdatam = csvread(filename + "-meta_x_t_times.csv");
names = ["50 ps"; "180 ps"]
plotmatched = 0;
%% Plot Raw Data
x = rawdatam(1, 1:189);
    
traw = rawdatam(2,:);
fig1 = figure(Name="Raw Data");
plot_D(rawdata, x, traw, 'Raw Data');
%% Remove Background
rawdata_avg = rawdata - mean(rawdata, 2);
fig2 = figure(Name="Background Removed")
titletext = "Background Removed";
hold on; colormap(gray(2^12));
imagesc(x,traw .*1e9,rawdata_avg);
title(titletext)
set(gca,'YDir','reverse');
ylabel('Time, ns');
xlabel('Scan x, m');
hold off
%% Low-Pass Filter
filteredData = lowpf(12e9, 50, rawdata_avg, 1/(x(2) - x(1)), 1/(traw(2) - traw(1)));
colormap(gray(2^12));
imagesc(x,traw .*1e9,filteredData);
set(gca,'YDir','reverse');
ylabel('Time, ns');
xlabel('Scan x, m');
%% Split into 2
data{1} = rawdata(1:end/2, :);
data{2} = rawdata(end/2:end, :);
t{1} = traw(1:end/2);
t{2} = traw(end/2:end);
%% Remove Background on each pulse Separately
for i = 1:height(names)
    data_avg{i} = data{i} - mean(data{i}, 2);
end
fig3 = figure(Name="Background Removal")o
tl = tiledlayout(2,1);
title(tl, "Background Removal")
for i = 1:height(names)
    nexttile;
    titletext = names(i);
    hold on; colormap(gray(2^12));
    imagesc(x,t{i} .*1e9,data_avg{i});
    title(titletext)
    set(gca,'YDir','reverse');
    ylabel('Time, ns');
    xlabel('Scan x, m');
    hold off
end
%% Low-Pass Filter on each pulse separately
kmax = [13, 13];
fmax = [4e9, 5e9];
for i = 1:height(names)
    data_avg_filtered{i} = lowpf(fmax(i), kmax(i), data_avg{i}, 1/(x(2) - x(1)), 1/(traw(2) - traw(1)));
end
fig4 = figure(Name="Low-Pass Filter")
tl = tiledlayout(2,1);
title(tl, "Low-Pass Filter")
for i = 1:height(names)
    nexttile;
    titletext = names(i);
    colormap(gray(2^12));
    imagesc(x,t{i} .*1e9,-data_avg_filtered{i});
    title(titletext)
    
    set(gca,'YDir','reverse');
    ylabel('Time, ns');
    xlabel('Scan x, m');
end
%% Matched filter
A_top = data_avg_filtered{1};
A_bottom = data_avg_filtered{2};
t_top = t{1}
t_bottom = t{1}
if plotmatched
    for degree = [1 2]
        
        data_MatchedFilter{1, degree} = matchedfilter(A_top, t_top, 50e-12, degree)
        data_MatchedFilter{2, degree} = matchedfilter(A_bottom, t_bottom, 180e-12, degree)
        
        figure(Name="Matched Filter " + degree);
        tl = tiledlayout(2,1);
        title(tl, "Matched Filter " + degree)
        
        for i = 1:height(names)
            nexttile;
            titletext = names(i);
            hold on; colormap(gray(2^10));
            imagesc(x,t{i}.*1e9,data_MatchedFilter{i});
            set(gca,'YDir','reverse');
            ylabel('Time, ns');
            xlabel('Scan x, m');
            title(titletext)
            hold off
        
        end
    
    end
end
%% Weiner Filter
ksize = 10;
useweiner = [1 1];
for i = 1:height(names)
    if useweiner(i)
        [data_Weiner{i}, nl] = wiener2(data_avg_filtered{i}, [ksize,ksize]); 
    else
        data_Weiner{i} = data_avg_filtered{i};
    end
end
fig5 = figure(Name="Weiner Filter ");
tl = tiledlayout(2,1);
title(tl, "Weiner Filter ")
for i = 1:height(names)
    nexttile;
    titletext = names(i);
    hold on; colormap(gray(2^10));
    imagesc(x,t{i}.*1e9,data_Weiner{i});
    set(gca,'YDir','reverse');
    ylabel('Time, ns');
    xlabel('Scan x, m');
    title(titletext)
end
%% Devolcolution weiner
freq = [60 180] * 1e-12;
SNR = 1.5e4;
for i = 1:height(names)
    psr = create_DOG_pulse(1, freq(i), t_top, t_top(end/2),1,1) % Define psr function as a differential pulse
    
    
    wavelet = psr.x
    wavelet = wavelet(end/2-200:end/2+200);
    deconv_data = data_avg_filtered{i};
    data_DevoncWiener{i} = deconvwnr(deconv_data, wavelet, SNR); % Apply Weiner filter to deconvolute the image
end
fig6 = figure(Name="Devonvolution Weiner Filter")
tl = tiledlayout(2,1);
title(tl, "Devonvolution Weiner Filter")
for i = 1:height(names)
    nexttile;
    titletext = names(i);
    hold on; colormap(gray(2^10));
    imagesc(x,t{i}.*1e9,-data_DevoncWiener{i});
    set(gca,'YDir','reverse');
    ylabel('Time, ns');
    xlabel('Scan x, m');
    title(titletext)
    hold off
end
%% recombine
data_DevoncWiener_combined = (data_DevoncWiener{1} + data_DevoncWiener{2}(1:end-1, :).*.2);
fig10 = figure(Name="Devonvolution Weiner Filter Combined")
title("Devonvolution Weiner Filter Combined")
hold on; colormap(gray(2^10));
imagesc(x,t{1}.*1e9,-data_DevoncWiener_combined);
ylim([0 10])
set(gca,'YDir','reverse');
ylabel('Time, ns');
xlabel('Scan x, m');
hold off
%%
ksize = 100;
useweiner = [1 1];
    [data_DevoncWiener_deblurr_combined, nl] = wiener2(data_DevoncWiener_combined, [ksize,ksize], .05e-11); 
% fig11 = figure(Name="Devonvolution Weiner Filter (Deblurred)");
title(tl, "Devonvolution Weiner Filter (Deblurred)")
    colormap(gray(2^10));
    imagesc(x,t{i}.*1e9,-data_DevoncWiener_deblurr_combined);
    set(gca,'YDir','reverse');
    ylabel('Time, ns');
    xlabel('Scan x, m');
%% Tyler recombine try
surface_adjust = 1;
if surface_adjust()
A_top = data_DevoncWiener{1};
A_bottom = data_DevoncWiener{2}(1:end-1, :);
for i = 1:size(A_top,1)
    for j = 1:size(A_top,2)
 %if i > length(t{1}) / 9
    A_top(i,j) = A_top(i,j) .* -exp(-i*0.003);
 %end
    end
end
    end
[top, bottom] = improvedAlignment(A_top,A_bottom,15);
tyler_combined  = (top.*5 + bottom) ./ 2;
fig10 = figure(Name="Devonvolution Weiner Filter Combined (Tyler)")
title("Devonvolution Weiner Filter Combined (Tyler)")
hold on; colormap(gray(2^10));
imagesc(x,t{1}.*1e9,-A_top);
ylim([0 10])
set(gca,'YDir','reverse');
ylabel('Time, ns');
xlabel('Scan x, m');
hold off
%% Optional de-blurring Weiner Filter
ksize = 40;
useweiner = [1 1];
for i = 1:height(names)
    if useweiner(i)
        [data_DevoncWiener_deblurr{i}, nl] = wiener2(data_DevoncWiener{i}, [ksize,ksize], 1e-12); 
    else
        data_DevoncWiener_deblurr{i} = data_DevoncWiener{i};
    end
end
fig7 = figure(Name="Devonvolution Weiner Filter (Deblurred)");
tl = tiledlayout(2,1);
title(tl, "Devonvolution Weiner Filter (Deblurred)")
for i = 1:height(names)
    nexttile;
    titletext = names(i);
    hold on; colormap(gray(2^10));
    imagesc(x,t{i}.*1e9,-data_DevoncWiener_deblurr{i});
    set(gca,'YDir','reverse');
    ylabel('Time, ns');
    xlabel('Scan x, m');
    title(titletext)
end
%% Optional Matched Filter After Weiner Filter
if plotmatched
    for degree = [1]
        k = 1;
        A_top = data_Weiner{1};
        A_bottom = data_Weiner{2};
        t_top = t{1}
        t_bottom = t{1}
        data_MatchedFilter{1, k} = matchedfilter(A_top, t_top, 50e-12, degree);
        data_MatchedFilter{2, k} = matchedfilter(A_bottom, t_bottom, 180e-12, degree);
        k  = k+1;
        
        fig8 = figure(Name="Weiner Filter + Matched Filter " + degree)
        tl = tiledlayout(2,1);
        title(tl, "Weiner Filter + Matched Filter " + degree)
        
        for i = 1:height(names)
            nexttile;
            titletext = names(i);
            hold on; colormap(gray(2^10));
            imagesc(x,t{i}.*1e9,data_MatchedFilter{i});
            set(gca,'YDir','reverse');
            ylabel('Time, ns');
            xlabel('Scan x, m');
            title(titletext)
        end
    
    end
end
%% Migration
migrate = 0;
if migrate
h=0;
er=3;
x_max = 1.5;
z_max = .25;
pixel_width = 1e-3;
xQ = 0:pixel_width:x_max;
zQ = pixel_width:pixel_width:z_max;
I_MF = BScan_kirchoff_migration4(data_avg_filtered{2}, h, er, t_top', x, -3e-2, 3e-2, xQ, zQ, 0);
DR=5; 
imgplot=20.*log10(abs(I_MF));
max1=max(imgplot(:));
figure(4); 
 hold on;
 colormap(jet);
imagesc(xQ,zQ,circshift(imgplot,0));
caxis([max1-DR max1]);
set(gca,'YDir','reverse');
xlim([min(xQ) max(xQ)]);
ylim([min(zQ) max(zQ)]);
ylabel('Depth, m');
xlabel('Distance, m');
xticks(0:0.25:1.5);
end
%%
function output = matchedfilter(A, tval, pulse_t, degree)
    MWFM = create_DOG_pulse(degree, pulse_t, tval, 2e-9,1,1);
    BKGR_MF=zeros(size(A)); %BKGR is background removed B-scan
    
    %perform 2d convolution (conv2 does not work... timing issue?)
    for k = 1:size(A,2)
        %flip the waveform due to conv() defi
        BKGR_MF(:,k)=conv(A(:,k),flip(MWFM.x),'same');
    end
    
    
    %throw away any negative correlation - not what are looking for
    BKGR_MF(BKGR_MF<0)=0;
    BKGR_MF = circshift(BKGR_MF,-1,1);
    output = BKGR_MF;
end
function plot_D(D,x,t,title_str)
    % create colorbar limits
    clim_abs_max = 250;
    threshold    = 0.99;
    D_sort = sort(D(:),'ascend');
    clims  = [D_sort(round(numel(D_sort)*(1-threshold))), D_sort(round(numel(D_sort)*threshold))];
    clims  = [max(-clim_abs_max, clims(1)),min(clim_abs_max, clims(2))];
    % figure();
    colormap(gray(2^12));
    imagesc(x, t.*1e9, D);
    colorbar;
    caxis(clims);
    title(title_str);
    xlabel("x (m)");
    ylabel("t (ns)");
end
function output = lowpf(f_max, kx_max, im, Kx, Fs)
B_fx = fftshift(fft(fftshift(fft(im,[],1),1),[],2),2);
kx = linspace(-Kx/2,Kx/2,size(im,2));
f = linspace(-Fs/2,Fs/2,size(im,1));
B_fx(abs(f)>f_max,:)=0;
B_fx(:,abs(kx)>kx_max)=0;
output = real(ifft(ifftshift(ifft(ifftshift(B_fx,2),[],2),1),[],1));
end
function [aligned_A, aligned_B] = improvedAlignment(A, B, max_shift)
    % Initialize variables
    [rows_A, cols_A] = size(A);
    [rows_B, cols_B] = size(B);
    max_corr = -Inf;
    best_shift = [0, 0];
    % Pad the B-scan matrices to allow for shifting
    padded_A = padarray(A, [max_shift, max_shift], 'replicate', 'both');
    padded_B = padarray(B, [max_shift, max_shift], 'replicate', 'both');
    % Iterate over possible shifts
    for row_shift = -max_shift:max_shift
        for col_shift = -max_shift:max_shift
            % Shift and compute cross-correlation
            shifted_B = circshift(padded_B, [row_shift, col_shift]);
            corr_val = sum(sum(padded_A .* shifted_B));
            % Update best shift if higher correlation is found
            if corr_val > max_corr
                max_corr = corr_val;
                best_shift = [row_shift, col_shift];
            end
        end
    end
    % Apply the best shift to align the B-scans
    aligned_B = circshift(B, best_shift);
    aligned_A = A;
end