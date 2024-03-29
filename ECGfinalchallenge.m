%% Calculate Breathing Rate & HR

%% Basic Signal Smoothing
% raw signals
ecg = data1(:,1); 
breath = data1(:,2);


% define baseline (mean) of each signal
ecg_baseline = mean(ecg);
breath_baseline = mean(breath);


% subtract baseline
ecg = ecg - ecg_baseline;
breath = breath - breath_baseline;

sample_rate = 2000; % Hz

%% Filtering
% bandpass filter for breathing signal
b_filt = designfilt("bandpassiir",FilterOrder=14, ...
    HalfPowerFrequency1=0.1,HalfPowerFrequency2=0.35, ...
    SampleRate=sample_rate);
breath_filtered = filtfilt(b_filt, breath);

% bandpass filter for ecg signal
ecg_filt = designfilt("bandpassiir",FilterOrder = 14, ...
    HalfPowerFrequency1 = 0.5,HalfPowerFrequency2 = 15, ...
    SampleRate = sample_rate);
ecg_filtered = filtfilt(ecg_filt, ecg);

%% Locating Breaths & Beats

% find each breath
[bpks,blocs] = findpeaks(breath_filtered);

% find each heartbeat (R wave) 
[hpks,hlocs] = findpeaks(ecg_filtered,MinPeakProminence=2);

%% Plotting

% create x axis for plotting
xaxis_points = zeros(length(breath_filtered),1);
for i = 1:length(breath_filtered)
    xaxis_points(i) = i;
end
xaxis_time = xaxis_points/sample_rate;

% plot breathing signal!
figure(1)
plot(xaxis_time,breath_filtered,'k')
hold on
scatter(blocs/2000,bpks,'m','filled')
title('Hardware & Software Filtered Breathing Signal',FontSize=15)
ylabel('Signal [mV]',FontSize=13)
xlabel('Time [s]',FontSize=13)
ax.FontSize = 16; 

% plot heart signal!
figure(2)
plot(xaxis_time,ecg_filtered,'k')
hold on
scatter(hlocs/2000,hpks,'m','filled')
title('Hardware & Software Filtered ECG Signal',FontSize=16)
ylabel('Signal [mV]',FontSize=14)
xlabel('Time [s]',FontSize=14)
ax.FontSize = 16; 

% plot heart signal!
figure(3)
plot(xaxis_time, ecg,'k')
hold on
title('Hardware Filtered ECG Signal',FontSize=16)
ylabel('Signal [mV]',FontSize=14)
xlabel('Time [s]',FontSize=14)
ax.FontSize = 16; 

% plot breathing signal!
figure(4)
plot(xaxis_time,breath,'k')
hold on
title('Hardware Filtered Breathing Signal', FontSize=16)
ylabel('Signal [mV]',FontSize=14)
xlabel('Time [s]',FontSize=14)
ax.FontSize = 16; 


