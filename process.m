%function process(y, fps)
clc
clear all
close all
pause(1)
load('data_test.mat');

% Parameters to play with
WINDOW_SECONDS = 6;             % [s] Sliding window length
BPM_SAMPLING_PERIOD = 0.5;      % [s] Time between heart rate estimations
BPM_L = 45; BPM_H = 200;        % [bpm] Valid heart rate range
FILTER_STABILIZATION_TIME = 1;  % [s] Filter startup transient
CUT_START_SECONDS = 0;          % [s] Initial signal period to cut off
FINE_TUNING_FREQ_INCREMENT = 1; % [bpm] Separation between test tones for smoothing
ANIMATION_SPEED_FACTOR = 2;     % [] This makes the animation run faster or slower than real time
figure(1);
plot(y);
title('Brightness Signal Computation');

% Build and apply input filter
[b, a] = butter(2, [(((BPM_L)/60)/fps*2) (((BPM_H)/60)/fps*2)]);
yf = filter(b, a, y);
y = yf((fps * max(FILTER_STABILIZATION_TIME, CUT_START_SECONDS))+1:size(yf, 2));
figure(2);
plot(y);
title('Filter Band-Pass');
% Some initializations and precalculations
num_window_samples = round(WINDOW_SECONDS * fps);
bpm_sampling_period_samples = round(BPM_SAMPLING_PERIOD * fps);
num_bpm_samples = floor((size(y, 2) - num_window_samples) / bpm_sampling_period_samples);
fcl = BPM_L / 60; fch = BPM_H / 60;
orig_y = y;
bpm = [];
bpm_smooth = [];

max_freq_plot_amplitude = 0;
max_time_plot_bpm = 100;
min_time_plot_bpm = 50;
MAX=0;

ROB = 18.5;
Hei=193;
Wei=98;
Agg=21;
Q=5;%COW
%Q = 4.5;%CEW
% Jurnal : 
% Parameter Estimation of the Arterial System
%Dubois D, Dubois EF. A formula to estimate the approximate surface area if height and weight be known. Arch Intern Med. 1916; 17:863-871.
BSA = 0.007184 * (power(Wei, 0.425)) * (power(Hei, 0.725));
sudah=0;
tsudah=0;
MAX_peak_v1=0;MAX_peak_v2=0;
for i=1:num_bpm_samples,
    
    % Isi sliding window dengan sinyal asli
    window_start = (i-1)*bpm_sampling_period_samples+1;
    ynw = orig_y(window_start:window_start+num_window_samples);
    % Gunakan jendela Hanning untuk membawa tepi ke nol. Dengan cara ini, tidak ada buatan
    % frekuensi tinggi muncul ketika sinyal diperlakukan sebagai periodik oleh
    % FFT
    y = ynw .* hann(size(ynw, 2))';
    if sudah==0
        sudah=1;
        figure(3)
        n=size(ynw, 2);
        X=[0:length(hann(size(ynw, 2)))-1]/30;
        Y=hann(size(ynw, 2))
        plot(X,Y);
        title('Hanning Windows');
        grid on;
        ylabel('Hanning Window Value');
        xlabel('Time (s)');
        pause(0.5);
    end
    gain = abs(fft(y));

    % Indeks FFT frekuensi di mana detak jantung manusia
    il = floor(fcl * (size(y, 2) / fps))+1; ih = ceil(fch * (size(y, 2) / fps))+1;
    index_range = il:ih;
    
    % Gambarkan amplitudo frekuensi yang diinginkan
    figure(4);
    subplot(2, 1, 1);
    hold off;
    
    fft_plot = plot((index_range-1) * (fps / size(y, 2)) * 60, gain(index_range), 'b', 'LineWidth', 2);
    hold on;    
    max_freq_plot_amplitude = max(max_freq_plot_amplitude, max(gain(index_range)));
    axis([BPM_L BPM_H 0 max_freq_plot_amplitude]);
    grid on;
    xlabel('Heart rate (BPM)');
    ylabel('Amplitude');
    title('Frequency analysis and time evolution of heart rate signal');

    % Temukan puncak dalam rentang frekuensi minat dan temukan yang tertinggi
    [pks, locs] = findpeaks(gain(index_range));
    
    [max_peak_v, max_peak_i] = max(pks) ;
    if (max_peak_v>MAX)
        MAX=max_peak_v;
    end
    max_f_index = index_range(locs(max_peak_i));
    bpm(i) = (max_f_index-1) * (fps / size(y, 2)) * 60;
    
    % Ratakan frekuensi puncak tertinggi dengan mencari frekuensi yang
    % "berkorelasi" terbaik dalam rentang resolusi di sekitar puncak
    freq_resolution = 1 / WINDOW_SECONDS;
    lowf = bpm(i) / 60 - 0.5 * freq_resolution;
    freq_inc = FINE_TUNING_FREQ_INCREMENT / 60;
    test_freqs = round(freq_resolution / freq_inc);
    Power = zeros(1, test_freqs);
    freqs = (0:test_freqs-1) * freq_inc + lowf;
    for h = 1:test_freqs,
        re = 0; im = 0;
        for j = 0:(size(y, 2) - 1),
            phi = 2 * pi * freqs(h) * (j / fps);
            re = re + y(j+1) * cos(phi);
            im = im + y(j+1) * sin(phi);
        end
        Power(h) = re * re + im * im;
    end
    [max_peak_v, max_peak_i] = max(Power);
    if max_peak_v>MAX_peak_v1
        MAX_peak_v1=max_peak_v;
        X1=(index_range-1) * (fps / size(y, 2)) * 60;
        Y1=gain(index_range);
        X2=freqs*60;
        Y2=sqrt(Power);
    end
    %max_peak_i
    bpm_smooth(i) = 60*freqs(max_peak_i);
    Beats=bpm(i);
    ET = (364.5 - 1.23 * Beats);
    % Jurnal : 
    % Parameter Estimation of the Arterial System
    SV = (-6.6 + (0.25 * (ET - 35)) - (0.62 * Beats) + (40.4 * BSA) - (0.51 * Agg));
    % Jurnal : 
    % Stroke Volume/Pulse Pressure Ratio and Cardiovascular
    % Risk in Arterial Hypertension
    PP = SV / ((0.013 * Wei - 0.007 * Agg - 0.004 * Beats) + 1.307);
    MPP = Q * ROB;
    % Plot amplitudo dalam interval fine-tuning
    hold on;
    smoothing_plot = plot(freqs*60, sqrt(Power), 'r', 'LineWidth', 2);

    SP(i) = (MPP + 3 / 2 * PP);%*2.0011;
    DP(i) = (MPP - PP / 3);%*9.1960;
    % Plot legend
    set(fft_plot, 'Displayname', 'FFT modulus');
    set(smoothing_plot, 'Displayname', 'Peak smoothing');
    legend('Location', 'NorthEast');
    
    % Plot BPM over time
    figure(4);
    subplot(2, 1, 2);
    t = (0:i-1) * ((size(orig_y, 2) / fps) / (num_bpm_samples - 1));
    hold off;
    plot(t, bpm_smooth, 'r', 'LineWidth', 2);
    max_time_plot_bpm = max(max_time_plot_bpm, max(bpm_smooth));
    min_time_plot_bpm = min(min_time_plot_bpm, min(bpm_smooth));
    axis([0 ((size(orig_y, 2)-1) / fps) min_time_plot_bpm max_time_plot_bpm]);
    grid on;
    xlabel('Time (s)');
    ylabel('Heart rate (BPM)');

    % Print speed factor over real time
    box = uicontrol('style', 'text');
    set(box, 'String', [num2str(ANIMATION_SPEED_FACTOR) 'x']);
    set(box, 'Position', [512-42, 7, 40, 38]);
    set(box, 'FontName', 'Ubuntu Condensed');
    set(box, 'ForegroundColor', hex2dec({'88' '88' '88'})/255);
    set(box, 'FontSize', 22);
    set(box, 'BackgroundColor', hex2dec({'cc' 'cc' 'cc'})/255);

    drawnow();  % Flush graphics    
    pause(BPM_SAMPLING_PERIOD / ANIMATION_SPEED_FACTOR);
        
end

disp(['Mean HR: ' num2str(mean(bpm_smooth)) ' bpm']);
 
A=floor(mean(SP(i)));
B=mean(SP(i)); 
%disp(['Mean DP: ' num2str((mean(DP(i))))]);
%disp(['Mean SP: ' num2str((mean(SP(i))))]);
if ((B-A)>0.49)
    disp(['Mean SP: ' num2str(ceil(mean(SP(i))))]);
else
    disp(['Mean SP: ' num2str(floor(mean(SP(i))))]);
end
A=floor(mean(DP(i)));
B=mean(DP(i)); 
%disp(['Mean DP: ' num2str((mean(DP(i))))]);
%disp(['Mean SP: ' num2str((mean(SP(i))))]);
if ((B-A)>0.49)
    disp(['Mean DP: ' num2str(ceil(mean(DP(i))))]);
else
    disp(['Mean DP: ' num2str(floor(mean(DP(i))))]);
end
[val,idx]=max(Y1);
X3=X1(idx-2:idx+1)
Y3=Y1(idx-2:idx+1)
figure(5)
plot(X1,Y1, 'b', 'LineWidth', 2);
hold on
plot(X2,Y2, 'r', 'LineWidth', 2);
title('Smoothing');
ylabel('Amplitude');
xlabel('Heart Rate (BPM)');
grid on
figure(6)
plot(X3,Y3, 'b', 'LineWidth', 2);
hold on
plot(X2,Y2, 'r', 'LineWidth', 2);
[val,idx]=max(Y2);
plot(X2(idx),val,'O')
text(X2(idx)+1,val,num2str(val))
title('Smoothing');
ylabel('Amplitude');
xlabel('Heart Rate (BPM)');
grid on
%end
