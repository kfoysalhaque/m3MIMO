% Clear all workspace variables except for sdr objects

% List of variables to keep
varsToKeep = {'sdr0', 'sdr1', 'tx_antennas_configs', 'config_idx'};

% Get the list of all variables in the workspace
allVars = whos;

% Construct the command to clear all variables except those in varsToKeep
clearCommand = 'clear ';
for k = 1:length(allVars)
    if ~ismember(allVars(k).name, varsToKeep)
        clearCommand = [clearCommand allVars(k).name ' '];
    end
end

% Execute the clear command
eval(clearCommand);

clear allVars clearCommand k

%% Generate the data

% Initialize a character for dynamic variable naming
character = '1';

sdrtx = sdr0;
sdrrx = sdr1;
analog_beamforming = true;   % all_channel_same_data

nFFT = 1024;	% number of FFT points
txPower = 8000;
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

constellation = [1+1j 1-1j -1+1j -1-1j]; 

txfd = zeros(nFFT, 8);
tx_antennas = [1, 2, 3, 4, 5, 6, 7, 8];

for i = 1:length(tx_antennas)
    antenna_index = tx_antennas(i);
    txfd(nFFT/2 + 1 + scMin:nFFT/2 + scMax, antenna_index) = tx_symbols;
end

txfd = fftshift(txfd, 1);
txtd = ifft(txfd, [], 1);

txtdMod = zeros(nFFT, sdrrx.nch);
for txIndex=1:sdrrx.nch
    txtdMod(:, txIndex) = txtd(:, txIndex) * exp(1j*txIndex*pi*sin(aod)); % Apply BF
end

txtd = sdrtx.applyCalTxArray(txtdMod);

% Normalize the energy of the tx array and scale with txPower.
txtd = txPower*txtd./max(abs(txtd));

nFFT = 1024;
nskip = 1024*3;	% skip ADC data
nbatch = 100;	% num of batches
sync_idx = 1;
channel_idx = 2;
data_idxs = 3;

%% Send the data to the DACs
sdrtx.send(txtd);

%% Collecting Data

% Initialize the collection variables
CFR_all = [];  % This will store all the collected CFRs

% Set up a figure to detect the space key press
hFig = figure;
set(hFig, 'KeyPressFcn', @(src, event) setappdata(hFig, 'flag', strcmp(event.Key, 'space')));

% Initialize the flag for stopping the loop
setappdata(hFig, 'flag', false);

disp('Press space to stop data collection');


while ~getappdata(hFig, 'flag')
    %% Receive data
    rxtd = sdrrx.recv(nFFT, nskip, nbatch);
    rxtd = sdrrx.applyCalRxArray(rxtd);

    %% Find the synchronization instant
    rxtd_sync = squeeze(rxtd(:, sync_idx, :));
    rxfd_sync = fft(rxtd_sync, [], 1);

    locs = zeros(length(tx_antennas), sdrrx.nch);  % timing offset
    for tx_idx = tx_antennas
        corrfd = txfd(:, tx_idx) .* conj(rxfd_sync);
        corrtd = ifft(corrfd, [], 1);
        [val, loc] = max(abs(corrtd), [], 1);
        loc = loc + 0; % add a timing offset to test
        locs(tx_idx, :) = loc;
    end

    % Remove outliers
    median_loc = round(median(locs, 'all'));
    locs(abs(locs-median_loc) > 5) = median_loc;
    loc_selected = round(mean(locs, 'all'));

    %% Estimate the channel through pilots
    rxtd_pilot = squeeze(rxtd(:, channel_idx, :));
    rxtd_pilot = [rxtd_pilot(nFFT-loc_selected+2:nFFT, :); ...
        rxtd_pilot(1:nFFT-loc_selected+1, :)];

    rxfd_pilot = fftshift(fft(rxtd_pilot, [], 1));
    txfd = fftshift(txfd, 1);

    CFR = zeros(nFFT, sdrrx.nch, sdrtx.nch);
    for rx_idx = 1:sdrrx.nch
        h_rx_idx = rxfd_pilot(:, rx_idx) ./ txfd;
        CFR(nFFT/2 + 1 + scMin:nFFT/2 + scMax, rx_idx, :) = ...
            h_rx_idx(nFFT/2 + 1 + scMin:nFFT/2 + scMax, :);
    end

    % Store the collected CFR
    CFR_all = cat(4, CFR_all, CFR);
end

disp('Data collection stopped.');

% Close the figure
close(hFig);

%%
name = ['gesture_CFR_' character '_' num2str(length(tx_antennas)) 'x' num2str(length(tx_antennas)) '.mat'];
save(name, "CFR_all", '-mat');

