% 2020.06.09
% Dianxin at UoE
% MMSE_Channel_Tap_Comb_Pilot is modified from the code of MIMO-OFDM
% WIRELESS COMMUNICATIONS WITH MATLAB

% Nfft is the size of IFFT/FFT
% h is the channel taps
% SNR is the signal to noise ratio of channel

% Because it is comb pilot, Pilot_location and Pilot_interval, and
% Received_Pilot, and Pilot_Value are required

function H_MMSE = MMSE_Channel_Tap_Comb_Pilot(Received_Pilot, Pilot_Value, Pilot_location, Nfft, Pilot_interval, SNR, h)

SNR_HEX = 10^(SNR / 10);
Np = size(Pilot_location, 2);
H_LS = Received_Pilot ./ Pilot_Value;

k = 0: length(h) - 1;
hh = h * h';
tmp = h .* conj(h) .* k;
r = sum(tmp) / hh;
r2 = tmp * k .'/hh;

tau_rms = sqrt(r2 - r^2);
df = 1/Nfft;
j2pi_tau_df = 1j * 2 * pi * tau_rms * df;
K1 = repmat([0 : Nfft - 1].', 1, Np);
K2 = repmat([0 : Np - 1], Nfft, 1);
rf = 1./(1 + j2pi_tau_df * (K1 - K2 * Pilot_interval));
K3 = repmat([0 : Np - 1].', 1, Np);
K4 = repmat([0 : Np - 1], Np, 1);
rf2 = 1./(1 + j2pi_tau_df * Pilot_interval * (K3 - K4));
Rhp = rf;
Rpp = rf2 + (eye(size(H_LS, 1)) / SNR_HEX);
H_MMSE = Rhp * pinv(Rpp) * H_LS;