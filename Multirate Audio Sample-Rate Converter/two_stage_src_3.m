Fs = 96e3;
t = 0:1/Fs:0.002;  % shorter time for plotting

% Example signal
f1 = 20e3;  
f2 = 30e3; 
x = sin(2*pi*f1*t) + sin(2*pi*f2*t);

%% Stage 1
L = 3;    % upsampling factor
M = 4;   % downsampling factor

% Step 1: Upsample (zero-stuff)
x_up = upsample(x, L);   

% Step 2: Design anti-imaging / anti-aliasing lowpass filter
Fs_up = Fs * L;
Fpass = 20e3;     
Fstop = 22e3;     
F_new_inter = Fs*L/2;

lpFilt = designfilt('lowpassfir', ...
    'PassbandFrequency', Fpass/(F_new_inter), ...
    'StopbandFrequency', Fstop/F_new_inter, ...
    'PassbandRipple', 0.01, ...
    'StopbandAttenuation', 100, ...
    'DesignMethod','kaiserwin');
b = lpFilt.Coefficients;
disp(['Stage 1 filter length = ', num2str(length(b))])

% Step 3: Filter
x_filt = filter(lpFilt, x_up);

% Step 4: Downsample
y = downsample(x_filt, M);

% Step 5: New sampling rate
Fs_new = Fs * (L/M);
fprintf('New Fs after Stage 1 = %.2f Hz\n', Fs_new);

%% Stage 2
L1 = 49;   % upsampling factor
M1 = 80;   % downsampling factor

% Step 1: Upsample
x_up_1 = upsample(y, L1);   

% Step 2: Design anti-imaging / anti-aliasing lowpass filter
F_new_inter = Fs_new*L1/2;
lpFilt1 = designfilt('lowpassfir', ...
    'PassbandFrequency', Fpass/(F_new_inter), ...
    'StopbandFrequency', Fstop/F_new_inter, ...
    'PassbandRipple', 0.01, ...
    'StopbandAttenuation', 100, ...
    'DesignMethod','kaiserwin');
b1 = lpFilt1.Coefficients;
disp(['Stage 2 filter length = ', num2str(length(b1))])

% Step 3: Filter
x_filt1 = filter(lpFilt1, x_up_1);

% Step 4: Downsample
y1 = downsample(x_filt1, M1);

% Step 5: New sampling rate
Fs_final = Fs_new * (L1/M1);
fprintf('New Fs after Stage 2 = %.2f Hz\n', Fs_final);

%% Plot
subplot(5,1,1);
stem(0:100, x(1:101), 'filled', 'MarkerSize', 3)
title('Original (96 kHz)')

subplot(5,1,2);
stem(0:443, x_up(1:444), 'filled', 'MarkerSize', 3)
title('Stage 1: After Upsampling (L=7)')

subplot(5,1,3);
stem(0:80, y(1:81), 'filled', 'MarkerSize', 3)
title('Stage 1: After Downsampling (M=16)')

subplot(5,1,4);
stem(0:100, x_up_1(1:101), 'filled', 'MarkerSize', 3)
title('Stage 2: After Upsampling (L=21)')


subplot(5,1,5);
stem(0:80, y1(1:81), 'filled', 'MarkerSize', 3)
title('Stage 2: After Downsampling (M=20)')
