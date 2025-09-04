% %% Add the folder containing +piradio to the MATLAB path.
addpath('../../');
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
nSymbols = 3;


nFFT = 1024;	% number of FFT points
txPower = 8000;
tx_antennas = [1, 2, 3, 4, 5, 6, 7, 8];

% scMin = 335;
% scMax = 350;
uncoded_bits = 228*4; % integer multiple of 4 (from the coding)
aod = 0;
aod = deg2rad(aod);

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

txfd_ant = zeros(nFFT, 8); 
% txfd = zeros(nFFT, nSymbols, 8); 

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
% Normalize the energy of the tx array and scale with txPower.
txtd = txPower * txtd ./ max(abs(txtd));

%% Send the data to the DACs for the first stage
sdrtx.send(txtd);

%% Receive data

% Configure the Switches
%sdr0.set_switches("normal");

nskip = nFFT * nSymbols * 3; % Skip ADC data
nbatch = 100; % Number of batches

for ind = 1:3
    rxtd = sdrrx.recv(nFFT * nSymbols, nskip, nbatch);
end

rxtd = sdrrx.applyCalRxArray(rxtd);
rxtd = rxtd(:, 1, :);
rxtd = reshape(rxtd, nFFT, nSymbols, 8);

%% Select packet indices for Synchronization and Channel Estimation for First Stage
sync_idx = 1;
channel_idx = 2;
data_idxs = 3;

%% Find the synchronization instant

rxtd_sync = squeeze(rxtd(:, sync_idx, :));
rxfd_sync = fftshift(fft(rxtd_sync, [], 1));

txfd_sync = fftshift(squeeze(txfd(:, sync_idx, :)), 1);

% Find the starting of the symbol
% it should be the same across antennas
% correlate each rx symbol (sum of all tx symbols) with each tx symbol

locs = zeros(size(tx_antennas, 2), sdrrx.nch);
for tx_idx = tx_antennas
    corrfd = txfd_sync(:, tx_idx) .* conj(rxfd_sync);
    corrtd = ifft(corrfd, [], 1);
    [val, loc] = max(abs(corrtd), [], 1);
    loc = loc + 0; % add a timing offset to test
    locs(tx_idx, :) = loc;
end

loc_selected = round(mean(locs, 'all'));

%% Estimate the channel through pilots
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

%% Equalize data through channel estimates
rxtd_data = squeeze(rxtd(:, data_idxs, :));
rxtd_data = [rxtd_data(nFFT-loc_selected+2:nFFT, :); ...
    rxtd_data(1:nFFT-loc_selected+1, :)];

rxfd_data = fftshift(fft(rxtd_data, [], 1));

t_hat = zeros(nFFT, size(tx_antennas, 2), sdrrx.nch);
for tx_idx = 1:size(tx_antennas, 2)
    h_idx = CFR(:, :, tx_idx);
    t_hat(:, tx_idx, :) = rxfd_data(:, :) ./ h_idx;
end

%% Plot
stream_idx = 1;
rx_idx = 8;
figure();
plot(real(t_hat(:, stream_idx, rx_idx)), ...
    imag(t_hat(:, stream_idx, rx_idx)), '.', 'MarkerSize', 20)
grid on;
ylim([-2 2]);
xlim([-2 2]);

%% BER
t_hat_ber = squeeze(t_hat(:, stream_idx, stream_idx));

rx_symbols =  t_hat_ber(nFFT/2 + 1 + scMin:nFFT/2 + scMax, 1);
[rx_bits, received_bits] = qpsk_demodulation(rx_symbols);
errors = sum(received_bits ~= input_bits.');
BER = errors / size(input_bits,1)


%% Precoding
% Multiplication for each of the subcarrier
CFR_permute = permute(CFR, [2, 3, 1]); % Nr-by-Nsts-by-Nst
[~,~,V] = pagesvd(CFR_permute, 'econ');  % CHECK
steeringMatrix = permute(V,[3 1 2]); % Permute to Nst-by-Ntx-by-Nsts
% Nst=K, Ntx=M, Nsts=Nss

% Zero-forcing precoding solution
for i = 1:nFFT
    % Channel inversion precoding
    h = squeeze(steeringMatrix(i,:,:));
    steeringMatrix(i,:,:) = h / (h'*h); 
end
% steeringMatrix should be K x M x Nss (Nss=M)


%% Prepare Second Stage Transmission with Precoding
nSymbols = 3;
txfd = zeros(nFFT, nSymbols, 8); % For 3 symbols in second stage

% Repeat txfd_ant along the second dimension
txfd = repmat(txfd_ant, [1, nSymbols, 1]);

% Reshape txfd to the desired shape 1024x3x8
txfd = reshape(txfd, [1024, nSymbols, 8]);


%% Apply precoding
txfd_precoded = zeros(nFFT, 1, 8);

txfd_i = squeeze(txfd(:, 3, :));
% for subcarrier = 1:nFFT
for subcarrier = (nFFT/2 + 1 + scMin):(nFFT/2 + scMax)
    txfd_precoded(subcarrier, 1, :) = squeeze(steeringMatrix(subcarrier, :, :)) * squeeze(txfd_i(subcarrier, :, :)).';
end


txfd(:,3,:) = txfd_precoded();
txfd(:,1,:) = txfd_precoded();

% % for subcarrier = 1:nFFT
% for subcarrier = 1:nFFT
%     txfd(subcarrier, 3, :) = squeeze(steeringMatrix(subcarrier, 1, :)) * squeeze(txfd(subcarrier, 3, :)).';
% end

txfd = fftshift(txfd, 1);
txtd = ifft(txfd, [], 1); % should be K x Nss (Nss=M)
txtd = reshape(txtd, nFFT * nSymbols, 8);

txtd = sdrtx.applyCalTxArray(txtd);
% Normalize the energy of the tx array and scale with txPower.
txtd = txPower * txtd ./ max(abs(txtd));

%% Send the data to the DACs for the second stage
sdrtx.send(txtd);

%% Receive data for the second stage
nskip = nFFT * nSymbols * 3;
nbatch = 100;

for ind = 1:3
    rxtd = sdrrx.recv(nFFT * nSymbols, nskip, nbatch);
end

rxtd = sdrrx.applyCalRxArray(rxtd);
rxtd = rxtd(:, 1, :);
rxtd = reshape(rxtd, nFFT, nSymbols, 8);

%% Synchronization for Second Stage
rxtd_sync = squeeze(rxtd(:, sync_idx, :));
rxfd_sync = fftshift(fft(rxtd_sync, [], 1));

txfd_sync = fftshift(squeeze(txfd(:, sync_idx, :)), 1);

% Find the starting of the symbol
% it should be the same across antennas
% correlate each rx symbol (sum of all tx symbols) with each tx symbol

locs = zeros(size(tx_antennas, 2), sdrrx.nch);
for tx_idx = tx_antennas
    corrfd = txfd_sync(:, tx_idx) .* conj(rxfd_sync);
    corrtd = ifft(corrfd, [], 1);
    [val, loc] = max(abs(corrtd), [], 1);
    loc = loc + 0; % add a timing offset to test
    locs(tx_idx, :) = loc;
end

loc_selected = round(mean(locs, 'all'));

%% Channel Estimation for Second Stage
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

%% Convert in the frequency domain
rxtd_data = squeeze(rxtd(:, data_idxs, :));
rxtd_data = [rxtd_data(nFFT-loc_selected+2:nFFT, :); ...
    rxtd_data(1:nFFT-loc_selected+1, :)];

rxfd_data = fftshift(fft(rxtd_data, [], 1));

%% Equalize data through channel estimates
t_hat = zeros(nFFT, size(tx_antennas, 2), sdrrx.nch);
for tx_idx = 1:size(tx_antennas, 2)
    h_idx = CFR(:, :, tx_idx);
    t_hat(:, tx_idx, :) = rxfd_data(:, :) ./ h_idx;
end
X_k = sum(t_hat, 3);

%% Plot
stream_idx = 1;
rx_idx = 1;
figure();
plot(real(t_hat(:, stream_idx, rx_idx)), ...
    imag(t_hat(:, stream_idx, rx_idx)), '.', 'MarkerSize', 20)
grid on;
ylim([-2 2]);
xlim([-2 2]);


%% BER
t_hat_bfm_NIC = t_hat(:, stream_idx, rx_idx);
t_hat_ber = squeeze(t_hat_bfm_NIC(:, 1));

rx_symbols =  t_hat_ber(nFFT/2 + 1 + scMin:nFFT/2 + scMax, 1);
[rx_bits, received_bits] = qpsk_demodulation(rx_symbols);
errors = sum(received_bits ~= input_bits.');
BER = errors / size(input_bits,1)

%% Equalize data through channel estimates with G
Nss_i = 8;
Nss = Nss_i;

I_Nss = eye(Nss_i, Nss);

% G matrix
G_k = zeros(nFFT, sdrtx.nch, sdrrx.nch);
% Final received symbol matrix
X_k = zeros(nFFT, sdrtx.nch);

% G_k for each subcarrier
for subcarrier = (nFFT/2 + 1 + scMin):(nFFT/2 + scMax)
    HW = squeeze(CFR(subcarrier, :, :)); % Combined H*W matrix for subcarrier
    Hw_C = HW';

    G = I_Nss * Hw_C / (HW * Hw_C + I_Nss);
    G_k(subcarrier,:,:) = G;

    Y = rxfd_data(subcarrier, :).';
    X_k(subcarrier, :) = Nss * G * Y;

end

%% Plot
stream_idx = 3;
figure();
plot(real(X_k(:, stream_idx)), ...
     imag(X_k(:, stream_idx)), '.', 'MarkerSize', 20);
grid on;
% ylim([-20000 20000]);
% xlim([-20000 20000]);


%% BER
t_hat_ber = X_k(:, stream_idx);

rx_symbols =  t_hat_ber(nFFT/2 + 1 + scMin:nFFT/2 + scMax, 1);
[rx_bits, received_bits] = qpsk_demodulation(rx_symbols);
errors = sum(received_bits ~= input_bits.');
BER = errors / size(input_bits,1)
%% Transmit processed data
sdrtx.send(zeros(nFFT, 8));

%% Receive data
sdrrx.recv(nFFT, nFFT * 3, 100);




