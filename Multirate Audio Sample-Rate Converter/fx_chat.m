%% ANC Headset DSP Core (Fx-LMS with Variable Step-Size)
% MATLAB R2023a-compatible, single-file script.
% Goal: 20 dB attenuation (100–1000 Hz) within 300 ms, <2% speech distortion (cepstral distance proxy).
% Notes:
% - Self-contained simulation (no audio files). Replace synthetic signals with real recordings as needed.
% - Uses feedforward ANC with one reference mic (x), one error mic at ear (e), and one loudspeaker.
% - Variable step-size NLMS based on normalized misalignment proxy and speech-protection (VAD-based).
% - Includes offline secondary-path identification (in sim we know S; in practice, replace with measured S_hat).

clear; close all; clc;

 %---------- 1) Parameters ----------
fs        = 8000;             % Sampling rate [Hz]; 8 kHz covers 100–1k Hz band
Ts        = 1/fs;
T         = 8.0;              % Total simulation [s]
N         = round(T*fs);

% Adaptive filter
Lw        = 128;              % ANC filter length (tune: 64–256)
mu_max    = 0.6;              % Max NLMS step-size (safe upper bound < 1)
mu_min    = 0.02;             % Floor step-size for tracking
leak      = 1e-4;             % Small leakage for robustness

% Secondary path (true, unknown to controller) and its estimate
Ls        = 64;               % Secondary-path length
sec_delay = 8;                % Samples delay in secondary path (models physical delay)

% Primary path (acoustic path of ambient noise to ear)
Lp        = 96;               % Primary path length

% Variable-step smoothing & power estimators
alpha_pwr = 0.90;             % EWMA for power (0.85–0.98)
alpha_vad = 0.95;             % EWMA for VAD envelope

% Performance measurement
band_lo   = 100;              % Hz
band_hi   = 1000;             % Hz
adapt_eval_start = 0.30;      % Evaluate attenuation after 300 ms

rng(2);

% ---------- 2) Synthetic signals ----------
% Reference noise (ambient before reaching ear)
% Road-traffic-like noise: low-passed, slightly 1/f-colored
w = randn(N,1);
[b_lp,a_lp] = butter(4, band_hi/(fs/2));        % LP to ~1kHz
w = filter(b_lp,a_lp,w);
w = filter([1 0.7],[1],w);                       % Add mild 1/f tilt
x = w / rms(w);                                  % Reference mic signal (normalized)

% Clean speech at ear (independent of reference). A simple voiced-like signal.
% You can replace with real speech vector `s_clean` sampled at fs.
t = (0:N-1)'/fs;
f0 = 140 + 10*sin(2*pi*0.5*t);                   % Slowly varying fundamental
s_clean = 0.2*sin(2*pi*cumsum(f0)/fs);
% Formants via resonators (rough approximation)
[bF1,aF1] = butter(2, [300 900]/(fs/2));
[bF2,aF2] = butter(2, [900 2500]/(fs/2));
s_clean = 0.6*filter(bF1,aF1,s_clean) + 0.4*filter(bF2,aF2,s_clean);
s_clean = 0.10*s_clean / max(abs(s_clean));       % Keep speech modest vs noise

% ---------- 3) Acoustic paths ----------
% Build primary path P(z) and secondary path S(z)
% Primary path: mildly resonant FIR (ear-cup leakage & reflections)
P = fir1(Lp-1, 0.35, kaiser(Lp,6));
P = P .* hanning(Lp)' ;
P = P / norm(P) * 0.8;                            % Gain scaling

% Secondary path S(z): delay + smooth LP FIR
S = zeros(Ls,1); S(sec_delay) = 1;                % Explicit delay
S = conv(S, fir1(Ls-1, 0.45));
S = S(1:Ls); S = S / norm(S) * 0.6;               % Normalize and scale

% Secondary-path estimate S_hat (with model error)
S_hat = S + 0.05*randn(size(S));                  % 5% model error
S_hat = S_hat / max(abs(S_hat));

% ---------- 4) Off-line secondary-path ID (optional demo) ----------
% In practice, you would identify S_hat by injecting a probe signal to the loudspeaker
% and measuring ear-mic response. Here we already set S_hat. You can uncomment to
% emulate an LMS identification step that produces S_hat.
%
% {  % Pseudo-code (commented):
% M_id = 1024; L_id = Ls; mu_id = 0.2; sid = zeros(L_id,1);
% u = randn(M_id,1); y = filter(S,1,u);  % loudspeaker drive and ear-mic
% Ubuf = zeros(L_id,1);
% for n=1:M_id
%   Ubuf = [u(n); Ubuf(1:end-1)];
%   e_id = y(n) - sid.'*Ubuf;
%   sid = sid + (mu_id/(1e-6+Ubuf.'*Ubuf))*e_id*Ubuf;
% end
% S_hat = sid;
% }

% ---------- 5) Fx-LMS ANC loop ----------
W = zeros(Lw,1);                % Adaptive filter coeffs (loudspeaker drive filter)
Xbuf = zeros(Lw,1);
Xfbuf = zeros(Lw,1);            % Filtered-x buffer

% Power trackers
Px = 1; Pe = 1; Py = 1; Pd = 1;   % EWMA powers

% Speech protection (simple VAD on error mic using band energy 200–3k Hz)
[b_vad,a_vad] = butter(2, [200, 3000]/(fs/2));
vad_env = 0; vad_thr = 1.5;       % Threshold in relative units

% Storage
y_hist  = zeros(N,1);            % Loudspeaker drive
es_hist = zeros(N,1);            % Secondary path output at ear from y
xpf_hist= zeros(N,1);            % Filtered-x instantaneous norm
mu_hist = zeros(N,1);
err     = zeros(N,1);            % Error mic signal (what the user hears)
noanc   = zeros(N,1);            % Error mic if ANC were off (for comparison)

for n = 1:N
    % --- Primary disturbance at ear (without ANC) ---
    % Noise through primary path
    xP = 0;
    for k=1:min(n,Lp)
        xP = xP + P(k)*x(n-k+1);
    end
    d_n = xP + s_clean(n);      % disturbance + desired speech

    % Reference buffer for adaptive filter
    Xbuf = [x(n); Xbuf(1:end-1)];

    % Anti-noise loudspeaker drive
    y_n = W.'*Xbuf;                     % y(n) = w^T x(n)
    y_hist(n) = y_n;

    % Secondary path to ear
    es_n = 0;
    for k=1:min(n,Ls)
        es_n = es_n + S(k)*y_hist(n-k+1);
    end
    es_hist(n) = es_n;

    % Error signal at ear mic (what the listener hears)
    e_n = d_n + es_n;              % sign convention: sum at ear; controller will learn cancellation
    err(n) = e_n;
    noanc(n) = d_n;                % baseline (ANC off)

    % Filtered-x signal using S_hat
    xfh = 0;
    for k=1:min(n,Ls)
        xfh = xfh + S_hat(k)*Xbuf(k);  % Using the same Xbuf ordering (approximation)
    end
    Xfbuf = [xfh; Xfbuf(1:end-1)];

    % --- Variable step-size (normalized misalignment proxy + VAD) ---
    % EWMA power updates
    Px = alpha_pwr*Px + (1-alpha_pwr)*(Xfbuf.'*Xfbuf);
    Pe = alpha_pwr*Pe + (1-alpha_pwr)*(e_n*e_n);

    % Output power proxy (anti-noise at ear)
    Py = alpha_pwr*Py + (1-alpha_pwr)*(es_n*es_n);
    Pd = alpha_pwr*Pd + (1-alpha_pwr)*(d_n*d_n);

    % Normalized misalignment proxy: higher when residual is large compared to anti-noise
    % rho in [0,1]; encourage larger mu when cancellation is poor, smaller mu when near convergence
    rho = min(1, max(0, (Pe) / (Pe + Py + 1e-12)));

    % VAD-based speech protection (reduce adaptation when speech dominates)
    vad_sample = filter(b_vad,a_vad,e_n);
    vad_env = alpha_vad*vad_env + (1-alpha_vad)*abs(vad_sample);
    speech_weight = 1 / (1 + (vad_env/vad_thr));   % ~1 when little speech; ~<1 when speech strong

    % NLMS step size with clamps
    mu = (mu_min + (mu_max-mu_min)*rho) * speech_weight;
    mu = max(mu_min, min(mu_max, mu));
    mu_hist(n) = mu;

    % --- Weight update (leaky NLMS) ---
    normXf = Px + 1e-8;           % Px already holds ||x_f||^2 EWMA
    W = (1 - leak)*W - (mu / normXf) * e_n * Xfbuf;

    % Debug: store x_f norm
    xpf_hist(n) = Px;
end

% ---------- 6) Metrics ----------
% Compute attenuation in band 100–1000 Hz after adaptation settles (>=300 ms)
idx_eval = round(adapt_eval_start*fs):N;
[ePxx,fgrid] = pwelch(noanc(idx_eval), hamming(2^12), [], [], fs);
[ePyy,~]     = pwelch(err(idx_eval),   hamming(2^12), [], [], fs);
band = (fgrid>=band_lo & fgrid<=band_hi);
att_dB = 10*log10( mean(ePxx(band)) / mean(ePyy(band)) );

% Estimate time-to-20dB: rolling bandpower ratio over time
win = round(0.20*fs);  % 200 ms window
step = round(0.02*fs);
ratios = [];
times  = [];
for i = 1:step:(N-win)
    seg = i:(i+win-1);
    [P0,fg] = pwelch(noanc(seg), hamming(512), [], [], fs);
    [P1,~]  = pwelch(err(seg),   hamming(512), [], [], fs);
    b = (fg>=band_lo & fg<=band_hi);
    ratios(end+1) = 10*log10(mean(P0(b))/mean(P1(b)));
    times(end+1)  = (i+win/2)/fs;
end
idx20 = find(ratios>20, 1, 'first');
if isempty(idx20)
    t20 = NaN;   % Did not reach 20 dB in the sim
else
    t20 = times(idx20);
end

% Cepstral distance-based speech distortion (proxy %)
% Compare cepstra of clean speech and the speech component at ear after ANC.
% We approximate "speech component at ear" by applying a spectral mask that
% deemphasizes sub-100 Hz and strong narrowband noise. This is a proxy – in
% lab use parallel path recording of speech alone.

% Frame-based MFCC-like cepstra (but we implement a simple real cepstrum)
frame_len = round(0.032*fs);        % 32 ms
hop       = round(0.010*fs);        % 10 ms

% High-pass to remove DC/very low-freq noise before cepstrum
[b_hp,a_hp] = butter(2, 80/(fs/2), 'high');
sp_ref = filter(b_hp,a_hp,s_clean);
sp_out = filter(b_hp,a_hp,err - (noanc - s_clean)); % heuristic isolation of speech

C_ref = real_cepstra(sp_ref, frame_len, hop);
C_out = real_cepstra(sp_out, frame_len, hop);

% Align frame counts
M = min(size(C_ref,2), size(C_out,2));
C_ref = C_ref(:,1:M); C_out = C_out(:,1:M);

% Cepstral distance per frame and normalized percent
cd = sqrt( sum( (C_ref - C_out).^2, 1 ) );
cd_pct = 100 * mean( cd ./ (sqrt(sum(C_ref.^2,1)) + 1e-12) );

% ---------- 7) Plots & Printouts ----------
figure; plot((0:N-1)/fs, [noanc err]); grid on;
xlabel('Time [s]'); ylabel('Amplitude'); legend('Ear signal (ANC off)','Ear signal (ANC on)');
title(sprintf('Time signals – Attenuation %.1f dB (100–1k Hz) after %.0f ms', att_dB, adapt_eval_start*1e3));

figure; plot(times, ratios); hold on; yline(20,'--'); grid on;
xlabel('Time [s]'); ylabel('Band attenuation [dB] (100–1k Hz)');
title('Adaptation trajectory');

figure; [P0,fg] = pwelch(noanc(idx_eval), hamming(2^12), [], [], fs);
[P1,~] = pwelch(err(idx_eval),   hamming(2^12), [], [], fs);
semilogx(fg,10*log10(P0),'-',fg,10*log10(P1),'-'); grid on;
xlabel('Frequency [Hz]'); ylabel('PSD [dB/Hz]'); xlim([50 2000]);
legend('ANC off','ANC on'); title('Post-adaptation PSD');

fprintf('\n=== PERFORMANCE SUMMARY ===\n');
fprintf('Mean attenuation (100–1000 Hz, after %.0f ms): %.1f dB\n', adapt_eval_start*1e3, att_dB);
if ~isnan(t20)
    fprintf('Time to reach 20 dB band attenuation: %.0f ms\n', t20*1e3);
else
    fprintf('Time to 20 dB band attenuation: NOT reached in this run\n');
end
fprintf('Cepstral-distance speech distortion (proxy): %.2f %% (target < 2%%)\n', cd_pct);

% ---------- 8) Tuning guidance ----------
% If attenuation < 20 dB or t20 > 300 ms:
% - Increase Lw (128 -> 192 or 256) and/or slightly raise mu_max (<=0.8) if stable.
% - Improve S_hat accuracy (lower modeling error, correct delay). A poor S_hat slows convergence.
% - Increase sec_delay to reflect true acoustic delay; mismatch hurts Fx-LMS.
% - Reduce alpha_pwr (faster power tracking) for quicker step-size response.
% If speech distortion > 2%%:
% - Lower vad_thr (more aggressive speech protection) or lower mu_max.
% - Increase alpha_vad (slower attack, more stable) to avoid adapting on speech.
% - Add a "freeze": if vad_env > vad_thr, set mu = mu_min.
% Stability tips:
% - Keep mu_max < 1 for NLMS; higher with precise S_hat but risk oscillation.
% - Leakage helps bound weights in non-stationary conditions.

% ---------- Helper: Real cepstra function ----------
function C = real_cepstra(x, frame_len, hop)
    % Returns low-quefrency real cepstral coefficients per frame (rows: quefrency, cols: frames)
    x = x(:);
    N = length(x);
    idx = 1;
    frames = {};
    while idx+frame_len-1 <= N
        seg = x(idx:idx+frame_len-1) .* hamming(frame_len);
        % Power spectrum
        X = fft(seg, 2^nextpow2(frame_len));
        L = log(abs(X).^2 + 1e-12);
        c = real(ifft(L));
        qmax = 20;  % keep 20 lowest quefrency coefficients (exclude c0)
        Cframe = c(2:(qmax+1));
        frames{end+1} = Cframe;
        idx = idx + hop;
    end
    if isempty(frames)
        C = zeros(20,0);
    else
        C = cell2mat(frames);
    end
end
