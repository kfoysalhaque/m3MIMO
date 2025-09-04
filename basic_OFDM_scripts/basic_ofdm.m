clear tx_antennas_configs config_idx tx_antennas

% List of tx_antennas configurations to test
tx_antennas_configs = {[1, 2, 3, 4, 5, 6, 7, 8], [1, 2, 3, 4]};

for config_idx = 1:length(tx_antennas_configs)
    tx_antennas = tx_antennas_configs{config_idx};

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

    rx_rotation = -60;

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

    % Normalize the energy of the tx array and scale with txPower.
    txtd = 20000*txtd./max(abs(txtd));


    txtdMod = zeros(nFFT, sdrrx.nch);
    for i = 1:length(tx_antennas)
        antenna_index = tx_antennas(i);
        txtdMod(:, antenna_index) = txtd;
    end

    txtdMod = sdrtx.applyCalTxArray(txtdMod);
    sdrtx.send(txtdMod);

    rxtd = sdrrx.recv(nread, nskip, ntimes);
    rxtd = sdrrx.applyCalRxArray(rxtd);

    %% Plot Receive Spectrum 

    naoa = 200;
    aoas = linspace(-1, 1, naoa);
    pArray = zeros(1, naoa);

    for iaoa = 1:naoa
        p = 0;
        aoa = aoas(iaoa);
        for itimes = 1:ntimes
            tdbf = zeros(nFFT, 1);
            for rxIndex=1:length(tx_antennas)
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
    ylim([-20 0])

    [maxValue, maxIndex] = max(pArray);
    disp(['Analog Beamforming Angle: ',  num2str(rad2deg(aoas(maxIndex)))]);


    %% SNR

    figure()
    rxtd_single_value = squeeze(rxtd(:, 1, 1:size(tx_antennas, 2)));
    a = sum(rxtd_single_value, 2);
    b = ((fftshift(abs(fft(a)))));
    plot(mag2db(b), 'LineWidth', 4);
    grid on; grid minor;
    startIndex = nFFT/2 + 1 + scMin;
    stopIndex = nFFT/2 + scMax;

    totSig2 = sum(abs(b(startIndex:stopIndex)) .* abs(b(startIndex:stopIndex)));
    S = totSig2 / (stopIndex-startIndex+1);
    totNoiA2 = sum(abs(b(1:startIndex-1)) .* abs(b(1:startIndex-1)));
    totNoiB2 = sum(abs(b(stopIndex+1:nFFT)) .* abs(b(stopIndex+1:nFFT)));
    totNoise2 = totNoiA2 + totNoiB2;
    N = totNoise2 / (2*startIndex - 3);
    snr = 10*log10(S/N);

    %% Save the Data
    rxtd_name = ['rxtd_' num2str(rx_rotation) '_' num2str(length(tx_antennas)) 'x' num2str(length(tx_antennas)) '.mat'];
    save(rxtd_name, "rxtd", '-mat');

    %% Stop transmitting and do a dummy read
    txtd = zeros(nFFT, sdrrx.nch);
    sdrtx.send(txtd);
    sdrrx.recv(nread, nskip, ntimes);
    pause(10);

end
