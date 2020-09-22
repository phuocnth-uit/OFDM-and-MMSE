% 2020.06.09
% Dianxin at UoE
% MMSE_Channel_Tap_Comb_Pilot is modified from the code of MIMO-OFDM
% WIRELESS COMMUNICATIONS WITH MATLAB

% Nfft is the size of IFFT/FFT
% h is the channel taps
% SNR is the signal to noise ratio of channel

% Frame_size is the thing that, how many blocks exsist within one coherent
% time. For example, if Frame_size is 2, then channel taps wont be changed
% over 2 blocks

% Pilot_location refers to which block is pilot and which block is data

% Received_Pilot and Pilot_Value are required

function H_MMSE_h = MMSE_Channel_Tap_Block_Pilot(Received_Pilot, Pilot_Value, Pilot_location, Nfft, Frame_size, SNR, h)

H_MMSE_h = zeros(Nfft, Frame_size);

SNR_HEX = 10^(SNR / 10);
Np = Nfft;
H_LS = Received_Pilot ./ Pilot_Value;
Nps = 1;

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
rf = 1./(1 + j2pi_tau_df * (K1 - K2 * Nps));
K3 = repmat([0 : Np - 1].', 1, Np);
K4 = repmat([0 : Np - 1], Np, 1);
rf2 = 1./(1 + j2pi_tau_df * Nps * (K3 - K4));
Rhp = rf;
Rpp = rf2 + (eye(length(H_LS)) / SNR_HEX);
H_MMSE = Rhp * pinv(Rpp) * H_LS;

for i_MMSE = 1 : size(H_MMSE, 1)
    H_MMSE_h(i_MMSE, :) = interp1(Pilot_location(:), H_MMSE(i_MMSE, :), 1:Frame_size, 'linear','extrap');
end