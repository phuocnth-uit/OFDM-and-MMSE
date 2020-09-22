% 2020.06.09
% Dianxin at UoE
% Ignore the Received_Signal_when_h_is_known

function [Unrecovered_Signal, Received_Signal_when_h_is_known] = OFDM_Receiver(Received_Signal, FFT_Size, Length_of_CP, Length_of_symbol, Received_Signal_when_h_is_known)

% S2P
Received_Signal = reshape(Received_Signal, Length_of_symbol, []);
Received_Signal_when_h_is_known = reshape(Received_Signal_when_h_is_known, Length_of_symbol, []); % When Channel State is known

% Remove CP
Received_signal_removed_CP = Received_Signal(Length_of_CP + 1 : end, :);
Received_signal_removed_CP_when_h_is_known = Received_Signal_when_h_is_known(Length_of_CP + 1 : end, :);

% FFT
Unrecovered_Signal = fft(Received_signal_removed_CP);
Unrecovered_Signal = (1 / sqrt(FFT_Size)) * Unrecovered_Signal; % check the power of Received_signal_removed_CP and Unrecovered_Signal, and the power is constant

Received_Signal_when_h_is_known = fft(Received_signal_removed_CP_when_h_is_known);
Received_Signal_when_h_is_known = (1 / sqrt(FFT_Size)) * Received_Signal_when_h_is_known;