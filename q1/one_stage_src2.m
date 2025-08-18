Fs = 96e3;
t = 0:1/Fs:0.002;  % shorter time for plotting

% Example signal
f1 = 20e3;  % tone 1
f2 = 30e3;  % tone 2
x = sin(2*pi*f1*t) + sin(2*pi*f2*t);

% Factors for sample rate conversion
L = 147;    % interpolation factor
M = 320;    % decimation factor

% Step 1: Upsample (zero-stuffing only, no filtering)
x_up = upsample(x, L);

% Step 2: Design anti-imaging / anti-aliasing filter
Fs_up = Fs * L;                  % New sampling rate after upsampling
Fpass = 20e3;                    % Passband cutoff (Hz)
Fstop = 22e3;                    % Stopband start (Hz)
F_nyq = Fs_up/2;                 % Nyquist after upsampling

lpFilt = designfilt('lowpassfir', ...
    'PassbandFrequency', Fpass/F_nyq, ...
    'StopbandFrequency', Fstop/F_nyq, ...
    'PassbandRipple', 0.01, ...
    'StopbandAttenuation', 100, ...
    'DesignMethod','kaiserwin');

fvtool(lpFilt)  % View filter characteristics
b = lpFilt.Coefficients;
disp(['Filter length: ', num2str(length(b))])

% Step 3: Filter (removes images from upsampling)
x_filt = filter(lpFilt, x_up);

% Step 4: Downsample
y = downsample(x_filt, M);

% Step 5: New sampling rate
Fs_new = Fs * (L/M);
fprintf('New Fs = %.2f Hz\n', Fs_new);

% Plot results
subplot(3,1,1);
stem(0:100, x(1:101), 'filled', 'MarkerSize', 3)
title('Original (96 kHz)')

subplot(3,1,2);
stem(0:443, x_up(1:444), 'filled', 'MarkerSize', 3)
title('After Upsampling (L=147, zero-stuffed)')

subplot(3,1,3);
stem(0:80, y(1:81), 'filled', 'MarkerSize', 3)
title(sprintf('After Downsampling (M=%d, Fs=%.1f kHz)', M, Fs_new/1e3))
