%% Lead Configurations
sampleRate = 2000;
x_axis1 = zeros(length(data1(:,1)),1);
for i = 1:length(data1)
    x_axis1(i) = i;
end
xaxis_time1 = x_axis1/sampleRate;

x_axis2 = zeros(length(data1(:,1)),1);
for i = 1:length(data1)
    x_axis2(i) = i;
end
xaxis_time2 = x_axis2/sampleRate;

x_axis3 = zeros(length(data1(:,1)),1);
for i = 1:length(data1)
    x_axis3(i) = i;
end
xaxis_time3 = x_axis3/sampleRate;

% Lead I
[filtered1_1,filtered2_1,filtered3_1] = FilterRawData(data1);
figure()
plot(x_axis1(20000:30000),filtered1_1(20000:30000))
hold on
title("Subject 1 Lead 1")
ylabel("Filtered ECG Signal [mV]")
xlabel("Time [s]")


% Lead II
[filtered1_2,filtered2_2,filtered3_2] = FilterRawData(data2);
figure()
plot(x_axis2(20000:30000),filtered1_2(20000:30000))
title("Subject 1 Lead 2")
ylabel("Filtered ECG Signal [mV]")
xlabel("Time [s]")


% Lead III
[filtered1_3,filtered2_3,filtered3_3] = FilterRawData(data3);
figure()
plot(x_axis3(20000:30000),filtered1_3(20000:30000))
ylim([-0.5 1.5])
title("Subject 1 Lead 3")
ylabel("Filtered ECG Signal [mV]")
xlabel("Time [s]")


