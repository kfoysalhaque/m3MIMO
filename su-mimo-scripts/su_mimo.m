% Specify the main folder
mainFolder = '/home/phd/MU-MIMO-Pi-Radio/new_radio';

% Generate the path string for the main folder and all its subfolders
allSubfolders = genpath(mainFolder);

% Add the paths to the MATLAB path
addpath(allSubfolders);

% Save the changes to the MATLAB path (optional)
savepath;

close all;

%% Configure the Switches
% sdr0.obsCtrl.configure(0);
% sdr1.obsCtrl.configure(0);


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

distance = '12m';
tx_antennas = [1, 2, 3, 4, 5, 6, 7, 8];
sdrtx = sdr0;
sdrrx = sdr1;
nSymbols = 3;
range = 1:8;



nFFT = 512;	% number of FFT points
txPower = 5000;
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

input_bits = randi([0, 1], 1, uncoded_bits)';

[tx_symbols_1, encoded_bits_1] = qpsk_modulation((randi([0, 1], 1, uncoded_bits)')); 
[tx_symbols_2, encoded_bits_2] = qpsk_modulation((randi([0, 1], 1, uncoded_bits)')); 
[tx_symbols_3, encoded_bits_3] = qpsk_modulation(input_bits); 


%% Save the Data
save_input = ['input_bits_1_' num2str(distance) '_' num2str(length(tx_antennas)) 'x' num2str(length(tx_antennas)) '.mat'];
save(save_input, "input_bits", '-mat');


%%
constellation = [1+1j 1-1j -1+1j -1-1j]; 


txfd_1 = zeros(nFFT, 8);
txfd_2 = zeros(nFFT, 8);
txfd_3 = zeros(nFFT, 8);


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

txfd = zeros(nFFT, nSymbols, 8);

txfd(:, 1, :) = txfd_1;
txfd(:, 2, :) = txfd_2;
txfd(:, 3, :) = txfd_3;

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


% 
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

%% Save the Data
CFR_unbeamformed = ['CFR_unbeamformed_' num2str(distance) '_' num2str(length(tx_antennas)) 'x' num2str(length(tx_antennas)) '.mat'];
save(CFR_unbeamformed, "CFR", '-mat');

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


%% Save the Data
equalized_t_hat = ['equalized_t_hat_' num2str(distance) '_' num2str(length(tx_antennas)) 'x' num2str(length(tx_antennas)) '.mat'];
save(equalized_t_hat, "t_hat", '-mat');


%% Plot
stream_idx = 1;
rx_idx = 8;
figure();
plot(real(t_hat(:, stream_idx, rx_idx)), ...
    imag(t_hat(:, stream_idx, rx_idx)), '.', 'MarkerSize', 20)
ylim([-2 2]);
xlim([-2 2]);
grid on;

%% BER

t_hat_no_bfm = mean(t_hat, 3);

t_hat_no_bfm = t_hat_no_bfm(:, range);

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

%% Precoding
% Multiplication for each of the subcarrier
CFR_permute = permute(CFR, [2, 3, 1]); % Nr-by-Nsts-by-Nst
[~,~,V] = pagesvd(CFR_permute, 'econ');
steeringMatrix = permute(V,[3 1 2]); % Permute to Nst-by-Ntx-by-Nsts
% Nst=K, Ntx=M, Nsts=Nss

% Zero-forcing precoding solution
for i = 1:nFFT
    % Channel inversion precoding
    h = squeeze(steeringMatrix(i,:,:));
    steeringMatrix(i,:,:) = (((h)*h')) / (h'*h);
    
end

% steeringMatrix should be K x M x Nss (Nss=M)

input_bits = randi([0, 1], 1, uncoded_bits)';

[tx_symbols_1, encoded_bits_1] = qpsk_modulation((randi([0, 1], 1, uncoded_bits)')); 
[tx_symbols_2, encoded_bits_2] = qpsk_modulation((randi([0, 1], 1, uncoded_bits)')); 
[tx_symbols_3, encoded_bits_3] = qpsk_modulation(input_bits); 

save_input = ['input_bits_2_' num2str(distance) '_' num2str(length(tx_antennas)) 'x' num2str(length(tx_antennas)) '.mat'];
save(save_input, "input_bits", '-mat');

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

txfd(:, 1, :) = txfd_1;
txfd(:, 2, :) = txfd_2;
txfd(:, 3, :) = txfd_3;

%% Apply precoding
txfd_precoded = zeros(nFFT, nSymbols, 8);
for symidx = 1:nSymbols
    txfd_i = squeeze(txfd(:, symidx, :));
    % for subcarrier = 1:nFFT
    for subcarrier = (nFFT/2 + 1 + scMin):(nFFT/2 + scMax)
        txfd_precoded(subcarrier, symidx, :) = squeeze(steeringMatrix(subcarrier, :, :)) * squeeze(txfd_i(subcarrier, :, :)).';
    end
end

txfd = fftshift(txfd_precoded, 1);
% txfd = fftshift(txfd, 1);
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


% 
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

%% Save the Data
CFR_beamformed = ['CFR_beamformed_' num2str(distance) '_' num2str(length(tx_antennas)) 'x' num2str(length(tx_antennas)) '.mat'];
save(CFR_beamformed, "CFR", '-mat');

%% Plot CFR
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
ylim([-2 2]);
xlim([-2 2]);
grid on;

%% BER

t_hat_no_bfm = mean(t_hat, 3);

t_hat_no_bfm = t_hat_no_bfm(:, range);
t_hat_a = mean(t_hat_no_bfm, 2);

rx_symbols =  t_hat_a(nFFT/2 + 1 + scMin:nFFT/2 + scMax, 1);
[rx_bits, received_bits] = qpsk_demodulation(rx_symbols);
errors = sum(received_bits ~= input_bits.');
BER_Precoding = errors / size(input_bits ,1)


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


%% Save the Data
beamformed_t_hat = ['beamformed_t_hat_' num2str(distance) '_' num2str(length(tx_antennas)) 'x' num2str(length(tx_antennas)) '.mat'];
save(beamformed_t_hat, "X_k", '-mat');


%% Plot
stream_idx = 3;
figure();
plot(real(X_k(:, stream_idx)), ...
     imag(X_k(:, stream_idx)), '.', 'MarkerSize', 20);
grid on;
ylim([-2 2]);
xlim([-2 2]);


%% BER
X_k = X_k(:, range);
X_k_mean = mean(X_k, 2);

rx_symbols =  X_k_mean(nFFT/2 + 1 + scMin:nFFT/2 + scMax, 1);
[rx_bits, received_bits] = qpsk_demodulation(rx_symbols);
errors = sum(received_bits ~= input_bits.');
BER_Final = errors / size(input_bits,1)


%% Equalized SNR
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

%% Transmit processed data
sdrtx.send(zeros(nFFT, 8));

%% Receive data
sdrrx.recv(nFFT, nFFT * 3, 100);

%% Specify the main folder
mainFolder = '/home/phd/MU-MIMO-Pi-Radio/new_radio';

% Generate the path string for the main folder and all its subfolders
allSubfolders = genpath(mainFolder);

% Add the paths to the MATLAB path
rmpath(allSubfolders);
% Save the changes to the MATLAB path (optional)
savepath;