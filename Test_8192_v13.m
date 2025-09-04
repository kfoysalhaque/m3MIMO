
% Specify the main folder
mainFolder = '/home/mentis/Github/MU-MIMO-Pi-Radio/new_radio';

% Generate the path string for the main folder and all its subfolders
allSubfolders = genpath(mainFolder);

% Add the paths to the MATLAB path
addpath(allSubfolders);

close all;

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
sdrtx = sdr0;
sdrrx = sdr1;
analog_beamforming = true;   % all_channel_same_data

nFFT = 512;	% number of FFT points
txPower = 8000;
% scMin = 335;
% scMax = 350;
uncoded_bits = 112*4; % integer multiple of 4 (from the coding)
aod = 0;
aod = deg2rad(aod);

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



nSymbols = 3; % NDP should only have two symbols, no data

txtd = zeros(nFFT, nSymbols, 8);



input_bits_3 = randi([0, 1], 1, uncoded_bits)';

[tx_symbols_1, encoded_bits_1] = qpsk_modulation((randi([0, 1], 1, uncoded_bits)')); 
[tx_symbols_2, encoded_bits_2] = qpsk_modulation((randi([0, 1], 1, uncoded_bits)')); 
[tx_symbols_3, encoded_bits_3] = qpsk_modulation(input_bits_3); 

constellation = [1+1j 1-1j -1+1j -1-1j]; 


txfd_1 = zeros(nFFT, 8);
txfd_2 = zeros(nFFT, 8);
txfd_3 = zeros(nFFT, 8);

tx_antennas = [1, 2, 3, 4, 5, 6, 7, 8];

for i = 1:length(tx_antennas)
    antenna_index = tx_antennas(i);
    txfd_1(nFFT/2 + 1 + scMin:nFFT/2 + scMax, antenna_index) = tx_symbols_1;
end


for i = 1:length(tx_antennas)
    antenna_index = tx_antennas(i);
    txfd_2(nFFT/2 + 1 + scMin:nFFT/2 + scMax, antenna_index) = tx_symbols_2;
end


for i = 1:length(tx_antennas)
    antenna_index = tx_antennas(i);
    txfd_3(nFFT/2 + 1 + scMin:nFFT/2 + scMax, antenna_index) = tx_symbols_3;
end


txfd_1 = fftshift(txfd_1, 1);
txfd_2 = fftshift(txfd_2, 1);
txfd_3 = fftshift(txfd_3, 1);


txtd_1 = ifft(txfd_1, [], 1);
txtd_2 = ifft(txfd_2, [], 1);
txtd_3 = ifft(txfd_3, [], 1);

% txtdMod = zeros(nFFT, sdrrx.nch);
% for txIndex=1:sdrrx.nch
%     txtdMod(:, txIndex) = txtd(:, txIndex) * exp(1j*txIndex*pi*sin(aod)); % Apply BF
% end

txtd_1 = sdrtx.applyCalTxArray(txtd_1);
txtd_2 = sdrtx.applyCalTxArray(txtd_2);
txtd_3 = sdrtx.applyCalTxArray(txtd_3);

txtd(:, 1, :) = txtd_1;
txtd(:, 2, :) = txtd_2;
txtd(:, 3, :) = txtd_3;

txtd_reshaped = reshape(txtd, [], 8);

txtd_sent = txPower*txtd_reshaped./max(abs(txtd_reshaped));

% Normalize the energy of the tx array and scale with txPower.
% txtd_1 = txPower*txtd_1./max(abs(txtd_1));
% txtd_2 = txPower*txtd_2./max(abs(txtd_2));
% txtd_3 = txPower*txtd_3./max(abs(txtd_3));


%% Send the data to the DACs
sdrtx.send(txtd_sent);

%% Receive data

% Configure the Switchestxfd_2
%sdr0.set_switches("normal");

nFFT =512;
nskip = nFFT * 4 * 3; % Skip ADC data
nbatch = 100; % Number of batches

%rxtd = sdrrx.recv(nFFT, nskip, nbatch);

for ind = 1:3
    rxtd = sdrrx.recv(nFFT * nSymbols, nskip, nbatch);
end

rxtd = sdrrx.applyCalRxArray(rxtd);
rxtd = rxtd(:, 1, :);
rxtd = reshape(rxtd, nFFT, nSymbols, 8);



%% Select packet indices for syncronization, channel estimation and data reception
sync_idx = 1;
channel_idx = 2;
data_idxs = 3;

%% Find the syncronization instant
rxtd_sync = squeeze(rxtd(:, sync_idx, :));
rxfd_sync = fft(rxtd_sync, [], 1);

% Find the starting of the symbol
% it should be the same across antennas
% correlate each rx symbol (sum of all tx symbols) with each tx symbol
locs = zeros(size(tx_antennas, 2), sdrrx.nch);  % timing offset
for tx_idx = tx_antennas
    corrfd = txfd_1(:, tx_idx) .* conj(rxfd_sync);
    corrtd = ifft(corrfd, [], 1);
    [val, loc] = max(abs(corrtd), [], 1);
    loc = loc+0; % add a timing offset to test
    locs(tx_idx, :) = loc;
end


% Remove outliers
median_loc = round(median(locs, 'all'));
locs(abs(locs-median_loc)>5) = median_loc;
loc_selected = round(mean(locs, 'all'));


%% Estimate the channel through pilots
rxtd_pilot = squeeze(rxtd(:, channel_idx, :));
% Compensate for the timing offset
rxtd_pilot = [rxtd_pilot(nFFT-loc_selected+2:nFFT, :); ...
    rxtd_pilot(1:nFFT-loc_selected+1, :)];

rxfd_pilot = fftshift(fft(rxtd_pilot, [], 1));
txfd_2 = fftshift(txfd_2, 1);
%%
figure(5)
for rx_idx = 1:sdrrx.nch
    subplot(4, 2, rx_idx);
    plot(abs(rxfd_pilot(:, rx_idx)), 'LineWidth', 2);
    title(strcat('FFT rx data ', num2str(tx_idx), '-', num2str(rx_idx)), 'FontSize', 15);
end
%%
CFR = zeros(nFFT, sdrrx.nch, sdrtx.nch);
for rx_idx = 1:sdrrx.nch
    h_rx_idx = rxfd_pilot(:, rx_idx) ./ txfd_2;
    CFR(nFFT/2 + 1 + scMin:nFFT/2 + scMax, rx_idx, :) = ...
        h_rx_idx(nFFT/2 + 1 + scMin:nFFT/2 + scMax, :);
end
%%
figure()
tx_idx = 1;
for rx_idx = 1:sdrrx.nch
    subplot(4, 2, rx_idx);
    plot(mag2db(abs(CFR(:, rx_idx, tx_idx))), 'LineWidth', 2);
    ylim([100 120]);
    xlim([0 nFFT]);
    title(strcat('CFR', num2str(tx_idx), '-', num2str(rx_idx)), 'FontSize', 20);
    set(gca, 'FontSize', 20); % Change 12 to your desired font size
end
%%
CIR = ifft(CFR, [], 1);
figure()
tx_idx = 1;
for rx_idx = 1:sdrrx.nch
    subplot(4, 2, rx_idx);
    stem(abs(CIR(:, rx_idx, tx_idx)), 'LineWidth', 2);
    title(strcat('CIR ', num2str(tx_idx), '-', num2str(rx_idx)), 'FontSize', 20);
    set(gca, 'FontSize', 20); % Change 12 to your desired font size
end

%% Equalize data through channel estimates
rxtd_data = squeeze(rxtd(:, data_idxs, :));
rxtd_data = [rxtd_data(nFFT-loc_selected+2:nFFT, :); ...
    rxtd_data(1:nFFT-loc_selected+1, :)];

rxfd_dataF = fftshift(fft(rxtd_data, [], 1));

t_hat = zeros(nFFT, size(tx_antennas, 2), sdrrx.nch);
for tx_idx = 1: sdrrx.nch
    h_idx = CFR(:, :, tx_idx);
    t_hat(:, tx_idx, :) = rxfd_dataF(:, :) ./ h_idx;
end
constellation = [1+1j 1-1j -1+1j -1-1j];

%% Plot

t_hat_no_bfm = mean(t_hat, 3);
stream_idx = 4;
rx_idx = 1;
figure();
plot(real(t_hat_no_bfm(:,  rx_idx)), ...
    imag(t_hat_no_bfm(:, rx_idx)), '.', 'MarkerSize', 20)
grid on;
ylim([-2 2]);
xlim([-2 2]);

%% BER

t_hat_no_bfm = mean(t_hat, 3);
t_hat_a = mean(t_hat, 2);

rx_symbols =  t_hat_no_bfm(nFFT/2 + 1 + scMin:nFFT/2 + scMax, 1);
[rx_bits, received_bits] = qpsk_demodulation(rx_symbols);
errors = sum(received_bits ~= input_bits_3.');
BER = errors / size(input_bits_3 ,1)

%% CFO compensation

%% Keeping only the data streams

 t_hat_data = t_hat(:, tx_antennas, : );

%% Combining
t_hat_comb = mean(t_hat_data, 3);
t_hat_all = mean(t_hat_comb, 2);

%% SNR

figure()
rxtd_single_value = squeeze(rxtd(:, 1, 1:size(tx_antennas, 2)));
a = sum(rxtd_single_value, 2);
b = ((fftshift(abs(fft(a)))));
plot(mag2db(b), 'LineWidth', 4);
grid on; grid minor;
startIndex = nFFT/2 + 1 + scMin;
stopIndex = nFFT/2 + 1 + scMax;

totSig2 = sum(abs(b(startIndex:stopIndex)) .* abs(b(startIndex:stopIndex)));
S = totSig2 / (stopIndex-startIndex+1);
totNoiA2 = sum(abs(b(1:startIndex-1)) .* abs(b(1:startIndex-1)));
totNoiB2 = sum(abs(b(stopIndex+1:nFFT)) .* abs(b(stopIndex+1:nFFT)));
totNoise2 = totNoiA2 + totNoiB2;
N = totNoise2 / (2*startIndex - 3);
snr = 10*log10(S/N)



%% Equalized SNR
n = 0;
sig = 0;
for sc = scMin:scMax-1
    if sc == 0
        continue;
    end

    % power domain calculations
    sig = sig + 2; % since it is QPSK
    x = abs(real(t_hat_all(nFFT/2 + 1 + sc))) - 1;
    y = abs(imag(t_hat_all(nFFT/2 + 1 + sc))) - 1;
    n = n + x^2 + y^2;
end
snr_equalized = 10*log10(sig/n)

%% BER
rx_symbols =  t_hat_all(nFFT/2 + 1 + scMin:nFFT/2 + scMax, 1);
[rx_bits, received_bits] = qpsk_demodulation(rx_symbols);
errors = sum(received_bits ~= input_bits_3.');
beamformed_BER = errors / size(input_bits_3,1)
figure();
plot(real(t_hat_all(:)), ...
    imag(t_hat_all(:)), '.', 'MarkerSize', 20)
grid on;
ylim([-2 2]);
xlim([-2 2]);
s = sprintf('ZF Equalized Symbols. SNR = %2.2f dB', snr_equalized);
title(s, 'FontSize', 20);
set(gca, 'FontSize', 20); % Change 12 to your desired font size

% Remove the paths from MATLAB's search path
rmpath(allSubfolders);




