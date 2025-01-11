clc;
close all;
clear all;
rotation = 60;
% aod = 0;
% aod = deg2rad(aod);
tx_antennas = [1, 2, 3, 4, 5, 6, 7, 8];
% tx_antennas = [1, 2, 3, 4];

addpath("E:\M-Flex\MU-MIMO-Pi-Radio\new_radio\beam_tracking");

analog_beamforming = true;   % all_channel_same_data

nFFT = 512;	% number of FFT points
txPower = 2000;
uncoded_bits = 112*4; % integer multiple of 4 (from the coding)


coding_rate = 7/4;
mod_rate = 1/2;
coded_bits_max = 400;
uncoded_bits_max = floor(coded_bits_max/mod_rate/coding_rate);

if uncoded_bits > uncoded_bits_max
    disp('ERROR, too many bits to transmit')
end

coded_bits = ceil(uncoded_bits*mod_rate*coding_rate);
scMin = -coded_bits/2;
scMax = coded_bits/2;

%% load the Data
unbeamformed_input = ['beam_tracking_input_bits_Nobfm_' num2str(rotation) '_' num2str(length(tx_antennas)) 'x' num2str(length(tx_antennas)) '.mat'];
load(unbeamformed_input, "input_bits");

unbeamformed_beam_tracking_t_hat = ['beam_tracking_t_hat_Nobfm_' num2str(rotation) '_' num2str(length(tx_antennas)) 'x' num2str(length(tx_antennas)) '.mat'];
load(unbeamformed_beam_tracking_t_hat, 't_hat');

t_hat_no_bfm = mean(t_hat, 3);
t_hat_a = mean(t_hat_no_bfm, 2);

%% Equalized SNR
n = 0;
sig = 0;
for sc = scMin:scMax-1
    if sc == 0
        continue;
    end

    % power domain calculations
    sig = sig + 2; % since it is QPSK
    x = abs(real(t_hat_a(nFFT/2 + 1 + sc))) - 1;
    y = abs(imag(t_hat_a(nFFT/2 + 1 + sc))) - 1;
    n = n + x^2 + y^2;
end
unbeamformed_snr = 10*log10(sig/n)

%% BER
rx_symbols =  t_hat_a(nFFT/2 + 1 + scMin:nFFT/2 + scMax, 1);
[rx_bits, received_bits] = qpsk_demodulation(rx_symbols);
errors = sum(received_bits ~= input_bits.');
unbeamformed_BER = errors / size(input_bits,1)
% figure();
% plot(real(t_hat_all(:)), ...
%     imag(t_hat_all(:)), '.', 'MarkerSize', 20)
% grid on;
% ylim([-2 2]);
% xlim([-2 2]);
% s = sprintf('ZF Equalized Symbols. SNR = %2.2f dB', unbeamformed_snr);
% title(s, 'FontSize', 20);
% set(gca, 'FontSize', 20); % Change 12 to your desired font size










%%
%%
%%
%%
%%
%%
%%
%%
%%   Beamformed Data


%% load the Data
beamformed_input = ['beam_tracking_input_bits_bfm_' num2str(rotation) '_' num2str(length(tx_antennas)) 'x' num2str(length(tx_antennas)) '.mat'];
load(beamformed_input, "input_bits");

beamformed_beam_tracking_t_hat = ['beam_tracking_t_hat_bfm_' num2str(rotation) '_' num2str(length(tx_antennas)) 'x' num2str(length(tx_antennas)) '.mat'];
load(beamformed_beam_tracking_t_hat, 't_hat');

t_hat_no_bfm = mean(t_hat, 3);
t_hat_a = mean(t_hat_no_bfm, 2);
%% Equalized SNR
n = 0;
sig = 0;
for sc = scMin:scMax-1
    if sc == 0
        continue;
    end

    % power domain calculations
    sig = sig + 2; % since it is QPSK
    x = abs(real(t_hat_a(nFFT/2 + 1 + sc))) - 1;
    y = abs(imag(t_hat_a(nFFT/2 + 1 + sc))) - 1;
    n = n + x^2 + y^2;
end
beamformed_snr = 10*log10(sig/n)

%% BER
rx_symbols =  t_hat_a(nFFT/2 + 1 + scMin:nFFT/2 + scMax, 1);
[rx_bits, received_bits] = qpsk_demodulation(rx_symbols);
errors = sum(received_bits ~= input_bits.');
beamformed_BER = errors / size(input_bits,1)


