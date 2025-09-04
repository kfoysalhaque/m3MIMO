% Add the folder containing +piradio to the MATLAB path.
addpath('../../');
close all;

% Clear all workspace variables except for sdr objects
varsToKeep = {'sdr0', 'sdr1'};
allVars = whos;
clearCommand = 'clear ';
for k = 1:length(allVars)
    if ~ismember(allVars(k).name, varsToKeep)
        clearCommand = [clearCommand allVars(k).name ' '];
    end
end
eval(clearCommand);
clear allVars clearCommand k

% Generate the data
sdrtx = sdr0;
sdrrx = sdr1;

nFFT = 1024;
nSymbols = 3;
txPower = 7000;
constellation = [1+1j 1-1j -1+1j -1-1j]; % QPSK constellation

txfd = zeros(nFFT, nSymbols, 8); % 3 symbols for the transmission
txfd_data = zeros(nFFT, 1, 8);
tx_antennas = [1, 2, 3, 4, 5, 6 , 7, 8];

uncoded_bits = 228*4; % integer multiple of 4 (from the coding)
coding_rate = 7/4;
mod_rate = 1/2;
coded_bits_max = 800;
uncoded_bits_max = floor(coded_bits_max/mod_rate/coding_rate);

if uncoded_bits > uncoded_bits_max
    disp('ERROR, too many bits to transmit')
end

coded_bits = ceil(uncoded_bits*mod_rate*coding_rate);
scMin = -coded_bits/2;
scMax = coded_bits/2;

input_bits = randi([0, 1], 1, uncoded_bits)';
[tx_symbols, encoded_bits] = qpsk_modulation(input_bits); 

% Modify only this part if you want to send identical symbols for all subcarriers
txfd_ant = zeros(nFFT, 8); 

for i = 1:length(tx_antennas)
    antenna_index = tx_antennas(i);
    txfd_ant(nFFT/2 + 1 + scMin:nFFT/2 + scMax, antenna_index) = tx_symbols;
end

% Repeat txfd_ant along the second dimension
txfd = repmat(txfd_ant, [1, nSymbols, 1]);

% Reshape txfd to the desired shape 1024x3x8
txfd = reshape(txfd, [1024, nSymbols, 8]);

txfd = fftshift(txfd, 1);
txtd = ifft(txfd, [], 1);
txtd = reshape(txtd, nFFT * nSymbols, 8);

txtd = sdrtx.applyCalTxArray(txtd);
txtd = txPower * txtd ./ max(abs(txtd));

sdrtx.send(txtd);

% Receive data
nskip = nFFT * nSymbols * 3; % Skip ADC data
nbatch = 100; % Number of batches

rxtd = sdrrx.recv(nFFT * nSymbols, nskip, nbatch);
rxtd = sdrrx.applyCalRxArray(rxtd);
rxtd = rxtd(:, 1, :);
rxtd = reshape(rxtd, nFFT, nSymbols, 8);

% Synchronization and Channel Estimation
sync_idx = 1;
channel_idx = 2;
data_idxs = 3;

rxtd_sync = squeeze(rxtd(:, sync_idx, :));
rxfd_sync = fftshift(fft(rxtd_sync, [], 1));

txfd_sync = fftshift(squeeze(txfd(:, sync_idx, :)), 1);

locs = zeros(size(tx_antennas, 2), sdrrx.nch);
for tx_idx = tx_antennas
    corrfd = txfd_sync(:, tx_idx) .* conj(rxfd_sync);
    corrtd = ifft(corrfd, [], 1);
    [val, loc] = max(abs(corrtd), [], 1);
    locs(tx_idx, :) = loc;
end

loc_selected = round(mean(locs, 'all'));
fprintf('Selected Synchronization Offset: %d\n', loc_selected);

% Estimate the channel through pilots
rxtd_pilot = squeeze(rxtd(:, channel_idx, :));
rxtd_pilot = [rxtd_pilot(nFFT - loc_selected + 2:nFFT, :); ...
    rxtd_pilot(1:nFFT - loc_selected + 1, :)];

rxfd_pilot = fftshift(fft(rxtd_pilot, [], 1));
txfd_pilot = fftshift(squeeze(txfd(:, channel_idx, :)), 1);

CFR = zeros(nFFT, sdrrx.nch, sdrtx.nch);
for rx_idx = 1:sdrrx.nch
    h_tmp = rxfd_pilot(:, rx_idx) ./ txfd_pilot;
    CFR(nFFT/2 + 1 + scMin:nFFT/2 + scMax, rx_idx, :) = ...
        h_tmp(nFFT/2 + 1 + scMin:nFFT/2 + scMax, :);
end

% Smoothing the CFR to reduce noise
CFR_smoothed = movmean(CFR, 5, 1);

% Plot Channel Frequency Response
figure;
subplot(2,1,1);
plot(abs(CFR_smoothed(:, 1, 1)));
title('Channel Frequency Response Magnitude (Smoothed)');
subplot(2,1,2);
plot(angle(CFR_smoothed(:, 1, 1)));
title('Channel Frequency Response Phase (Smoothed)');

% Equalize data through channel estimates
rxtd_data = squeeze(rxtd(:, data_idxs, :));
rxtd_data = [rxtd_data(nFFT-loc_selected+2:nFFT, :); ...
    rxtd_data(1:nFFT - loc_selected + 1, :)];

rxfd_data = fftshift(fft(rxtd_data, [], 1));

t_hat = zeros(nFFT, size(tx_antennas, 2), sdrrx.nch);
for tx_idx = 1:size(tx_antennas, 2)
    h_idx = CFR_smoothed(:, :, tx_idx);
    t_hat(:, tx_idx, :) = rxfd_data(:, :) ./ h_idx;
end

% Plot Equalized Symbols
t_hat_no_bfm = mean(t_hat, 3);
t_hat_ber = squeeze(t_hat_no_bfm(:, 1));

rx_symbols =  t_hat_ber(nFFT/2 + 1 + scMin:nFFT/2 + scMax, 1);
figure;
plot(real(rx_symbols), imag(rx_symbols), '.', 'MarkerSize', 20);
title('Equalized Symbols');
grid on;
ylim([-2 2]);
xlim([-2 2]);

% BER Calculation
[rx_bits, received_bits] = qpsk_demodulation(rx_symbols);
errors = sum(received_bits ~= input_bits.');
BER = errors / size(input_bits,1);
fprintf('Bit Error Rate (BER): %f\n', BER);

% Modulation function
function [tx_symbols, encoded_bits] = qpsk_modulation(bits)
    encoded_bits = encode(bits, 7, 4, 'hamming/binary');
    symbols_I = 2*encoded_bits(1:2:end) - 1; % Map bits to I-channel: 0 -> -1, 1 -> 1
    symbols_Q = 2*encoded_bits(2:2:end+1) - 1; % Map bits to Q-channel: 0 -> -1, 1 -> 1
    tx_symbols = symbols_I + 1j * symbols_Q;
end

% Demodulation function
function [rx_bits, corrected_bits] = qpsk_demodulation(rx_symbols)
    rx_bits_I = real(rx_symbols) > 0; % Demodulate I-channel: 1 -> 1, 0 -> 0
    rx_bits_Q = imag(rx_symbols) > 0; % Demodulate Q-channel: 1 -> 1, 0 -> 0
    rx_bits = zeros(1, 2 * length(rx_symbols));
    rx_bits(1:2:end) = rx_bits_I;
    rx_bits(2:2:end) = rx_bits_Q;
    corrected_bits = decode(rx_bits, 7, 4, 'hamming/binary');
end
