    % In this demo, we assume that sdrrx and sdrtx open, and is fully
% calibrated. Look at the calibration demo to make sure this is done. In
% the minimum, the timing and phase offsets need to be calibrated.

% Transmit a wideband signal from one channel on the TX. On the RX, capture
% samples, and apply the calibrations. Then, apply BF vectors for a set of
% AoA values. Plot them out.
%% Clear all workspace variables except for sdr objects

% List of variables to keep
varsToKeep = {'sdr0', 'sdr1', 'tx_antennas_configs', 'config_idx', 'tx_antennas'};

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

%%
sdrtx = sdr0;
sdrrx = sdr1;
rotation = '60';
tx_antennas = [1, 2, 3, 4, 5, 6, 7, 8];

nFFT = 1024;
nread = nFFT;
nskip = nFFT*3;
ntimes = 100;
txfd = zeros(nFFT, 1);
constellation = [1+1j 1-1j -1+1j -1-1j];

sdrrx.set_switches('off');

scMin = -400; scMax = 400;

for scIndex = scMin:scMax
    txfd(nFFT/2 + 1 + scIndex) = constellation(randi(4));
end

txfd = fftshift(txfd);
txtd = ifft(txfd);
m = max(abs(txtd));
txtd = txtd / m * 20000;

refTxIndex = 1;
txtdMod = zeros(nFFT, sdrrx.nch);
txtdMod(:, refTxIndex) = txtd;

txtdMod = sdrtx.applyCalTxArray(txtdMod);
sdrtx.send(txtdMod);

rxtd = sdrrx.recv(nread, nskip, ntimes);
rxtd = sdrrx.applyCalRxArray(rxtd);

%% Save the Data
rx_bfm = ['rx_bfm_' num2str(rotation) '_' num2str(length(tx_antennas)) 'x' num2str(length(tx_antennas)) '.mat'];
save(rx_bfm, "rxtd", '-mat');

naoa = 101;
aoas = linspace(-1.57, 1.57, naoa);
pArray = zeros(1, naoa);

for iaoa = 1:naoa
    p = 0;
    aoa = aoas(iaoa);
    for itimes = 1:ntimes
        tdbf = zeros(nFFT, 1);
        for rxIndex=1:sdrrx.nch
            td = rxtd(:,itimes,rxIndex);
            tdbf = tdbf + td * exp(1j*rxIndex*pi*sin(aoa)); % Apply BF Vec
        end % rxIndex
        fd = fftshift(fft(tdbf));
        p = p + sum(abs(fd( nFFT/2 + 1 + scMin : nFFT/2 + 1 + scMax)));
    end %itimes
    pArray(iaoa) = p;
end % iaoa

% Plot
pArray = pArray / max(pArray);
figure(3); clf;
plot(rad2deg(aoas), mag2db(pArray));
xlabel('Angle of Arrival (Deg)');
ylabel('Power (dB)');
grid on; grid minor;
ylim([-15 0])

[maxValue, maxIndex] = max(pArray);
disp(['Analog Beamforming Angle: ',  num2str(rad2deg(aoas(maxIndex)))]);


% Stop transmitting and do a dummy read
txtd = zeros(nFFT, sdrrx.nch);
sdrtx.send(txtd);
sdrrx.recv(nread, nskip, ntimes);

% Clear workspace variables
clear aoa aoas fd iaoa naoa p pArray refTxIndex td tdbf txtdMod;
clear ans itimes m nFFT nread nskip ntimes rxIndex rxtd;
clear scIndex txfd txtd constellation scMax scMin;
