% Test script for calcularMetricas.m
clc; clear; close all;

% 1. Create Synthetic Data
fs = 1e6;
N = 1024; % Freq bins
M = 100;  % Time steps (segments)
f = linspace(1.41e9, 1.43e9, N); % Freq vector

% Noise (Gaussian in Voltage -> Rayleigh in Magnitude -> Exponential in Power)
noise_power_linear = 1e-9 * ones(N, M); % -90 dBm base
% Add variance (Rayleigh fading / Chi-square 2 dof for power)
% For simplicity, let's just make Log-Normal or similar variance in dB
% Proper way: Complex Gaussian Noise
noise_real = randn(N, M);
noise_imag = randn(N, M);
noise_linear = (noise_real.^2 + noise_imag.^2) * 1e-9; % Mean power ~ 2e-9

% Add Signal (Continuous Wave at center)
idx_sig = N/2;
signal_power_linear = 100 * 1e-9; % 100x noise power (+20 dB)
% Add signal to the noise at that bin
noise_linear(idx_sig, :) = noise_linear(idx_sig, :) + signal_power_linear;

% Convert to dBm
% 1 mW = 0 dBm. 1e-9 W = 1e-6 mW = -60 dBm? No.
% 1e-9 W = -60 dBm?
% P(dBm) = 10*log10(P(W) / 1mW).
% If P(W) = 1e-9. P(mW) = 1e-6. 10*log10(1e-6) = -60 dBm.
P_dBm_matrix = 10*log10(noise_linear * 1000);

% Expected:
% Noise Floor ~ -60 + 10*log10(2) approx -57 dBm (due to Chi2 mean being 2*var/2... wait)
% Let's just check the values.
% Signal ~ 100*Noise. So Signal should be ~20 dB above noise.

% 2. Run calcularMetricas
metrics = calcularMetricas(f, P_dBm_matrix, 1:M, []);

disp('--- Test Results ---');
disp(['Calculated Difference: ' num2str(metrics.Diferencia_Potencia_dBm) ' dB']);
disp(['Peak Power: ' num2str(metrics.P_max_dBm) ' dBm']);
disp(['Noise Est: ' num2str(metrics.P_ruido_dBm) ' dBm']);

% 3. Check "Visual" difference (Average Profile)
P_profile_avg = mean(P_dBm_matrix, 2);
[p_peak_avg, ~] = max(P_profile_avg);
p_noise_avg = median(P_profile_avg);
diff_avg = p_peak_avg - p_noise_avg;

disp(['Average Profile Difference: ' num2str(diff_avg) ' dB']);

% 4. Check "Max Hold" difference (what the script does)
P_profile_max = max(P_dBm_matrix, [], 2);
[p_peak_max, ~] = max(P_profile_max);
p_noise_max = median(P_profile_max);
diff_max = p_peak_max - p_noise_max;

disp(['Max Profile Difference: ' num2str(diff_max) ' dB']);
