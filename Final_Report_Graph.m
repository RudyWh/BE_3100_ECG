%% Day 2 A1

%% Passive High-pass Filter

% define constants
R = 3.183e6; % ohms
C = 1000e-12; % F

% transfer function
    % https://en.wikipedia.org/wiki/High-pass_filter
numerator1 = [R*C,0];   % RCs
denominator1 = [R*C,1]; % RCs+1
sys1 = tf(numerator1,denominator1);


[mag1, phase1, wout1] = bode(sys1);
wout1 = wout1/(2*pi);  % convert to Hz

mag_plot1 = zeros(size(mag1,3),1);
for k = 1:size(mag1,3)
    mag_plot1(k,1) = mag1(:,:,k);
end

phase_plot1 = zeros(size(phase1,3),1);
for k = 1:size(phase1,3)
    phase_plot1(k,1) = phase1(:,:,k);
end

% experimental data
passive_data = readmatrix("Passive_Bode_Fixed.xlsx");
freq_p = passive_data(2:21,1);
mag_p = passive_data(2:21,2);
phase_p = passive_data(2:21,6);


% BODE plot
% subplot(2,1,1);
subplot(1,2,1);
semilogx(wout1(1:34),mag_plot1(1:34)); % in dB
hold on
xline(50)
hold on;
scatter(freq_p,mag_p);
hold off;
title("Magnitude response");
xlabel("log(Frequency) (log(Hz))");
ylabel("Magnitude (dB)");
legend("Model","","Experimental Data");

subplot(1,2,2);
semilogx(wout1(1:34),phase_plot1(1:34));
hold on
xline(50)
hold on;
scatter(freq_p,phase_p);
hold off;
title("Phase response");
xlabel("log(Frequency) (log(Hz))");
ylabel("Phase (deg)");
legend("Model","","Experimental Data");

sgtitle("Passive high-pass filter");


%% Active Low-pass Filter

% define constants
R1 = 50e3; % ohms
R2 = 50e3; % ohms
R4 = 50e3; % ohms
RF1 = 3.3e6; % ohms
RF2 = 3.3e6; % ohms
RG = 50e3; % ohms
RQ = 44.2e3; % ohms
C1 = 1000e-12; % F
C2 = 1000e-12; % F

Omega_n_sqd = R2/(R1*RF1*RF2*C1*C2);  
Q = (1+(R4/RQ)) * (1/((1/R1)+(1/R2)+(1/RG))) * sqrt((RF1*C1)/(R1*R2*RF2*C2));
A_LP = R1/RG;

% transfer function
numerator2 = [0,0,A_LP*Omega_n_sqd];
denominator2 = [1,sqrt(Omega_n_sqd)/Q,Omega_n_sqd];
sys2 = tf(numerator2,denominator2);

[mag2, phase2, wout2] = bode(sys2);
wout2 = wout2/(2*pi);
mag_plot2 = zeros(size(mag2,3),1);
for k = 1:size(mag2,3)
    mag_plot2(k,1) = mag2(:,:,k);
end
phase_plot2 = zeros(size(phase2,3),1);
for k = 1:size(phase2,3)
    phase_plot2(k,1) = phase2(:,:,k);
end

% experimental data
active_data = readmatrix("Active_Bode_Fixed.xlsx");
freq_a = active_data(:,1);
mag_a = active_data(:,2);
phase_a = active_data(:,6);

% BODE plot
% subplot(2,1,2);
figure;
subplot(1,2,1);
semilogx(wout2(1:34),mag_plot2(1:34));
hold on
xline(50)
hold on;
scatter(freq_a,mag_a);
hold off;
title("Magnitude response");
xlabel("log(Frequency) (log(Hz))");
ylabel("Magnitude (dB)");
legend("Model","","Experimental Data");

subplot(1,2,2);
semilogx(wout2(1:34),phase_plot2(1:34));
hold on
xline(50)
hold on;
scatter(freq_a,phase_a);
hold off;
title("Phase response");
xlabel("log(Frequency) (log(Hz))");
ylabel("Phase (deg)");
legend("Model","","Experimental Data");

sgtitle("Active low-pass filter");


