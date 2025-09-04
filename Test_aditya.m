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
txPower = 7000;

nFFT = 1024;	% number of FFT points
txPower = 20000;
tx_antennas = [1, 2, 3, 4, 5, 6, 7, 8];
nSymbols = 3;

txfd = zeros(nFFT, nSymbols, 8);

f1 = 50;
f2 = 100;
f3 = 150;
f4 = 200;
txfd(nFFT/2+f1+1, 1, 1) = 1;
%txfd(nFFT/2-f1+1, 1, 1) = 1;
txfd(nFFT/2+f2+1, 2, 1) = 1;
%txfd(nFFT/2-f2+1, 2, 1) = 1;
txfd(nFFT/2+f3+1, 3, 1) = 1;
%txfd(nFFT/2-f3+1, 3, 1) = 1;
% txfd(nFFT/2+f4+1, 4, 1) = 1;

figure(); plot(real(txfd(:, 1, 1))); grid on; xlim([1, 1024]); title("TX SYMB 1")
figure(); plot(real(txfd(:, 2, 1))); grid on; xlim([1, 1024]); title("TX SYMB 1")
% figure(); plot(real(txfd(:, 3, 1))); grid on; xlim([1, 1024]); title("TX SYMB 1")

txfd = fftshift(txfd, 1);
txtd = (ifft(txfd, [], 1)); % should be K x Nss (Nss=M)
% figure(); plot(real(txtd(:, 1, 1)));
% figure(); plot(real(txtd(:, 2, 1)));
% figure(); plot(real(txtd(:, 3, 1)));

txtd_resh = reshape(txtd, nFFT * nSymbols, 8);

txtd_resh = txPower * txtd_resh ./ max(abs(txtd_resh));

%% Send the data to the DACs for the second stage
sdrtx.send(txtd_resh);

%% Receive data for the second stage
nskip = nFFT * nSymbols * 3;
nbatch = 100;

for ind = 1:3
    rxtd = sdrrx.recv(nFFT * nSymbols, nskip, nbatch);
end

rxtd = sdrrx.applyCalRxArray(rxtd);
rxtd = squeeze(rxtd(:, 1, :));
rxtd = reshape(rxtd, nFFT, nSymbols, 8);

% figure(); plot(real(rxtd(:, 1, 1)));
% figure(); plot(real(rxtd(:, 2, 1)));
% figure(); plot(real(rxtd(:, 3, 1)));

rxfd = fftshift(fft(rxtd, [], 1));

figure(); plot(abs(rxfd(:, 1, 1))); grid on; xlim([1, 1024]); title("RX SYMB 1")
figure(); plot(abs(rxfd(:, 2, 1))); grid on; xlim([1, 1024]); title("RX SYMB 2")
figure(); plot(abs(rxfd(:, 3, 1))); grid on; xlim([1, 1024]); title("RX SYMB 3")
% figure(); plot(abs(rxfd(:, 4, 1))); grid on; xlim([1, 1024]); title("RX SYMB 4")
