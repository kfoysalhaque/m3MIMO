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

for ns = 1:2
    txfd(nFFT/2 + 1 + scMin:nFFT/2 + scMax, ns, tx_antennas) = ...
        constellation(randi(size(constellation, 2), scMax-scMin, size(tx_antennas, 2)));
end

for i = 1:length(tx_antennas)
    antenna_index = tx_antennas(i);
    txfd_data(nFFT/2 + 1 + scMin:nFFT/2 + scMax, antenna_index) = tx_symbols;
end

txfd(:, 3, :) = txfd_data;

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

% Synchronize each symbol independently
loc_selected = zeros(nSymbols, 1);

for sym_idx = 1:nSymbols
    rxtd_sync = squeeze(rxtd(:, sym_idx, :));
    rxfd_sync = fftshift(fft(rxtd_sync, [], 1));
    txfd_sync = fftshift(squeeze(txfd(:, sym_idx, :)), 1);
    
    locs = zeros(size(tx_antennas, 2), sdrrx.nch);
    for tx_idx = tx_antennas
        corrfd = txfd_sync(:, tx_idx) .* conj(rxfd_sync);
        corrtd = ifft(corrfd, [], 1);
        [val, loc] = max(abs(corrtd), [], 1);
        locs(tx_idx, :) = loc;
    end
    
    loc_selected(sym_idx) = round(mean(locs, 'all'));
    fprintf('Selected Synchronization Offset for symbol %d: %d\n', sym_idx, loc_selected(sym_idx));
end

% Apply synchronization to the entire received data based on the selected offset for each symbol
for sym_idx = 1:nSymbols
    rxtd(:, sym_idx, :) = circshift(rxtd(:, sym_idx, :), -loc_selected(sym_idx) + 1, 1);
end

% Channel Estimation using pilots from the second symbol
rxtd_pilot = squeeze(rxtd(:, 2, :));
rxfd_pilot = fftshift(fft(rxtd_pilot, [], 1));
txfd_pilot = fftshift(squeeze(txfd(:, 2, :)), 1);

CFR = zeros(nFFT, sdrrx.nch, sdrtx.nch);
for rx_idx = 1:sdrrx.nch
    h_tmp = rxfd_pilot(:, rx_idx) ./ txfd_pilot;
    CFR(nFFT/2 + 1 + scMin:nFFT/2 + scMax, rx_idx, :) = ...
        h_tmp(nFFT/2 + 1 + scMin:nFFT/2 + scMax, :);
end

% Equalize data through channel estimates using data from the third symbol
rxtd_data = squeeze(rxtd(:, 3, :));
rxfd_data = fftshift(fft(rxtd_data, [], 1));

t_hat = zeros(nFFT, size(tx_antennas, 2), sdrrx.nch);
for tx_idx = 1:size(tx_antennas, 2)
    h_idx = CFR(:, :, tx_idx);
    t_hat(:, tx_idx, :) = rxfd_data(:, :) ./ h_idx;
end

% Plot
stream_idx = 1;
rx_idx = 8;
figure();
plot(real(t_hat(:, stream_idx, rx_idx)), ...
    imag(t_hat(:, stream_idx, rx_idx)), '.', 'MarkerSize', 20)
grid on;
ylim([-2 2]);
xlim([-2 2]);

% BER Calculation
t_hat_no_bfm = mean(t_hat, 3);
t_hat_ber = squeeze(t_hat_no_bfm(:, 1));

rx_symbols = t_hat_ber(nFFT/2 + 1 + scMin:nFFT/2 + scMax, 1);
[rx_bits, received_bits] = qpsk_demodulation(rx_symbols);
errors = sum(received_bits ~= input_bits.');
BER = errors / size(input_bits, 1);
fprintf('Bit Error Rate (BER): %f\n', BER);

% Debugging: Print input bits and received bits
disp('Input Bits:');
disp(input_bits(1:100)); % Print first 100 input bits for brevity

disp('Received Bits:');
disp(received_bits(1:100)); % Print first 100 received bits for brevity

% Modulation function
function [tx_symbols, encoded_bits] = qpsk_modulation(bits)
    encoded_bits = encode(bits, 7, 4, 'hamming/binary');
    symbols_I = 2 * encoded_bits(1:2:end) - 1; % Map bits to I-channel: 0 -> -1, 1 -> 1
    symbols_Q = 2 * encoded_bits(2:2:end + 1) - 1; % Map bits to Q-channel: 0 -> -1, 1 -> 1
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
