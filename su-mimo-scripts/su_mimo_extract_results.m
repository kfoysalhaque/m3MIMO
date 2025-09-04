%% Clear all workspace variables except for sdr objects

% List of variables to keep
varsToKeep = {'sdr0', 'sdr1'};

% Get the list of all variables in the workspace
allVars = whos;

% Construct the command to clear all variables except those in varsToKeep
clearCommand = 'clear ';
for k = 1:length(allVars)
    if ~ismember(allVars(k).name, varsToKeep)
        clearCommand = [clearCommand allVars(k).name ' '];
    end
end

% Execute the clear commandtxfd_ant
eval(clearCommand);

clear allVars clearCommand k

%% Generate the data
clc;
close all;
clear all;
addpath('E:\M-Flex\MU-MIMO-Pi-Radio\new_radio\su-mimo_data')
distance = '1m';
tx_antennas = [1, 2, 3, 4, 5, 6, 7, 8];
range = 1:8;

nFFT = 512;	% number of FFT points
txPower = 8000;
uncoded_bits = 112*4;
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




save_input = ['input_bits_1_' num2str(distance) '_' num2str(length(tx_antennas)) 'x' num2str(length(tx_antennas)) '.mat'];
load(save_input, 'input_bits');


%% Equalized BER

equalized_t_hat = ['su-mimo_data/equalized_t_hat_' num2str(distance) '_' num2str(length(tx_antennas)) 'x' num2str(length(tx_antennas)) '.mat'];
load(equalized_t_hat, 't_hat');

%% Equalized CFR
CFR_unbeamformed = ['su-mimo_data/CFR_unbeamformed_' num2str(distance) '_' num2str(length(tx_antennas)) 'x' num2str(length(tx_antennas)) '.mat'];
load(CFR_unbeamformed, 'CFR');

%% Plot
% % stream_idx = 1;
% % rx_idx = 8;
% t_hat_no_bfm = mean(t_hat, 3);
% t_hat_a = mean(t_hat_no_bfm, 2);
% figure();
% plot(real(t_hat_a(:, 1)), ...
%     imag(t_hat_a(:, 1)), '.', 'MarkerSize', 20)
% ylim([-2 2]);
% xlim([-2 2]);
% grid on;
% set(gca, 'FontSize', 20);
% title('Received QPSK constellations',  'FontSize', 20);
% 
% %% Save plot
% plot_filename = ['su-mimo_nobfm_Constellation_' num2str(distance) '_' num2str(length(tx_antennas)) 'x' num2str(length(tx_antennas))];
% print(gcf, [plot_filename '.eps'], '-depsc');

%% BER

t_hat_no_bfm = mean(t_hat, 3);
t_hat_a = mean(t_hat_no_bfm, 2);

rx_symbols =  t_hat_a(nFFT/2 + 1 + scMin:nFFT/2 + scMax, 1);
[rx_bits, received_bits] = qpsk_demodulation(rx_symbols);
errors = sum(received_bits ~= input_bits.');
BER_Initial = errors / size(input_bits ,1)

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
snr_equalized_no_Beamform = 10*log10(sig/n)

%% Plot CFR
% figure()
% tx_idx = 1;
% for rx_idx = 1:8
%     subplot(4, 2, rx_idx);
%     plot(mag2db(abs(CFR(:, rx_idx, tx_idx))), 'LineWidth', 2);
%     ylim([100 120]);
%     xlim([0 nFFT]);
%     title(strcat('CFR', num2str(tx_idx), '-', num2str(rx_idx)), 'FontSize', 20);
%     set(gca, 'FontSize', 20); % Change 12 to your desired font size
% end
% 
% %% Save plot
% plot_filename = ['su-mimo_nobfm_CFR_' num2str(distance) '_' num2str(length(tx_antennas)) 'x' num2str(length(tx_antennas))];
% print(gcf, [plot_filename '.eps'], '-depsc');



%%
%%
%%
%%
%%
%%
%% Beamformed BER

save_input = ['su-mimo_data/input_bits_2_' num2str(distance) '_' num2str(length(tx_antennas)) 'x' num2str(length(tx_antennas)) '.mat'];
load(save_input, 'input_bits');

beamformed_t_hat = ['su-mimo_data/beamformed_t_hat_' num2str(distance) '_' num2str(length(tx_antennas)) 'x' num2str(length(tx_antennas)) '.mat'];
load(beamformed_t_hat, "X_k");

%% Beamformed CFR
CFR_beamformed = ['su-mimo_data/CFR_beamformed_' num2str(distance) '_' num2str(length(tx_antennas)) 'x' num2str(length(tx_antennas)) '.mat'];
load(CFR_beamformed, 'CFR');

%% Plot
% X_k_a = mean(X_k, 2);
% figure();
% plot(real(X_k_a(:, 1)), ...
%      imag(X_k_a(:, 1)), '.', 'MarkerSize', 20);
% grid on;
% ylim([-2 2]);
% xlim([-2 2]);
% set(gca, 'FontSize', 20);
% title('Received QPSK constellations',  'FontSize', 20);
% %% Save plot
% plot_filename = ['su-mimo_bfm_Constellation_' num2str(distance) '_' num2str(length(tx_antennas)) 'x' num2str(length(tx_antennas))];
% print(gcf, [plot_filename '.eps'], '-depsc');

%% BER
X_k = X_k(:, range);
X_k_mean = mean(X_k, 2);

rx_symbols =  X_k_mean(nFFT/2 + 1 + scMin:nFFT/2 + scMax, 1);
[rx_bits, received_bits] = qpsk_demodulation(rx_symbols);
errors = sum(received_bits ~= input_bits.');
BER_Final = errors / size(input_bits,1)


%% Beamformed SNR
n = 0;
sig = 0;
for sc = scMin:scMax-1
    if sc == 0
        continue;
    end

    % power domain calculations
    sig = sig + 2; % since it is QPSK
    x = abs(real(X_k_mean(nFFT/2 + 1 + sc))) - 1;
    y = abs(imag(X_k_mean(nFFT/2 + 1 + sc))) - 1;
    n = n + x^2 + y^2;
end
snr_equalized_Beamformed = 10*log10(sig/n)


%% Plot CFR
figure('Position', [100, 100, 900, 450]);

tx_idx = 1;
for rx_idx = 1:4
    subplot(2, 2, rx_idx);
    plot(mag2db(abs(CFR(:, rx_idx, tx_idx))), 'LineWidth', 2);
    ylim([100 120]);
    yticks([100 110 120]);
    yticklabels({'100', '110', '120'});
    % xlim([0 nFFT]);
    xlim([(nFFT/2 + 1 + scMin-5) (nFFT/2 + scMax-5)]);
    title(strcat('TX stream-', num2str(tx_idx), '-', 'RX stream-',  num2str(rx_idx)), 'FontSize', 16);
    set(gca, 'FontSize', 16); % Change 12 to your desired font size
    xlabel('Sub-channel', FontSize=16);
    ylabel('Magnitude', FontSize=16);
    grid on;
end

% Save plots
plot_filename = ['su-mimo_bfm_CFR_' num2str(distance) '_' num2str(length(tx_antennas)) 'x' num2str(length(tx_antennas))];
print(gcf, [plot_filename '.eps'], '-depsc');
