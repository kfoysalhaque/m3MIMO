close all;
clc;
clear all;

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

% Execute the clear command
eval(clearCommand);

clear allVars clearCommand k
close all

FntSize = 25;
rx_rotations = [0, 45, -45];
configs = {[1, 2, 3, 4, 5, 6, 7, 8]};

% Define common x-ticks and y-ticks for directional response pattern
common_xticks_dir = -90:30:90;
common_yticks_dir = -20:5:0;

% Define common y-ticks for receive power spectrum
common_yticks_pw = -50:10:0;

% Define figure size (width, height) in pixels
figure_width = 600;
figure_height = 500;

for i = 1:length(rx_rotations)
    rx_rotation = rx_rotations(i);
    for j = 1:length(configs)
        tx_antennas = configs{j};

        %% Variables
        nFFT = 1024;
        nread = nFFT;
        nskip = nFFT*3;
        ntimes = 100;
        txfd = zeros(nFFT, 1);
        constellation = [1+1j 1-1j -1+1j -1-1j];
        scMin = -400; scMax = 400;
        
        %% Load rxtd
        rxtd_name = ['basic_OFDM_data/rxtd_' num2str(rx_rotation) '_' num2str(length(tx_antennas)) 'x' num2str(length(tx_antennas)) '.mat'];
        load(rxtd_name, 'rxtd');
        
        %% SNR
        rxtd_single_value = squeeze(rxtd(:, 1, 1:size(tx_antennas, 2)));
        a = sum(rxtd_single_value, 2);
        b = ((fftshift(abs(fft(a)))));
        grid on; grid minor;
        startIndex = nFFT/2 + 1 + scMin;
        stopIndex = nFFT/2 + 1 + scMax;
        
        totSig2 = sum(abs(b(startIndex:stopIndex)) .* abs(b(startIndex:stopIndex)));
        S = totSig2 / (stopIndex-startIndex+1);
        totNoiA2 = sum(abs(b(1:startIndex-1)) .* abs(b(1:startIndex-1)));
        totNoiB2 = sum(abs(b(stopIndex+1:nFFT)) .* abs(b(stopIndex+1:nFFT)));
        totNoise2 = totNoiA2 + totNoiB2;
        N = totNoise2 / (2*startIndex - 3);
        snr = 10*log10(S/N);
        
        %% Plot Directional Response Pattern
        naoa = 200;
        aoas = linspace(1.5708, -1.5708, naoa);
        pArray = zeros(1, naoa);
        
        for iaoa = 1:naoa
            p = 0;
            aoa = aoas(iaoa);
            for itimes = 1:ntimes
                tdbf = zeros(nFFT, 1);
                for rxIndex = 1:length(tx_antennas)
                    td = rxtd(:, itimes, rxIndex);
                    tdbf = tdbf + td * exp(1j * rxIndex * pi * sin(aoa)); % Apply BF Vec
                end % rxIndex
                fd = fftshift(fft(tdbf));
                p = p + sum(abs(fd(nFFT/2 + 1 + scMin : nFFT/2 + 1 + scMax)));
            end % itimes
            pArray(iaoa) = p;
        end % iaoa
        
        %% Plot
        pArray = pArray / max(pArray);
        figure('Position', [100, 100, figure_width, figure_height]); % Set figure size
        clf;
        plot(rad2deg(aoas), mag2db(pArray), 'LineWidth', 5); % Increase LineWidth
        xlabel('Angle of Arrival (Deg)', 'FontSize', FntSize);
        ylabel('Power (dB)', 'FontSize', FntSize);
        set(gca, 'FontSize', FntSize, 'XTick', common_xticks_dir, 'YTick', common_yticks_dir);
        grid on; grid minor;
        set(gca, 'LineWidth', 0.5); % Increase grid line width
        ylim([-20 0]);


        %% Add SNR as an annotation in the figure
        str = sprintf('SNR = %.2f dB', snr);

        if rx_rotation > 0
            annotation('textbox', [0.2, 0.8, 0.1, 0.1], 'String', str, 'FitBoxToText', 'on', 'BackgroundColor', 'white', 'FontSize', FntSize);
        else
            annotation('textbox', [0.5, 0.8, 0.1, 0.1], 'String', str, 'FitBoxToText', 'on', 'BackgroundColor', 'white', 'FontSize', FntSize);
        end

        plot_filename = ['dir_res_pattern_' num2str(rx_rotation) '_' num2str(length(tx_antennas)) 'x' num2str(length(tx_antennas))];
        print(gcf, [plot_filename '.eps'], '-depsc');

        [maxValue, maxIndex] = max(pArray);
        disp(['Analog Beamforming Angle: ', num2str(rad2deg(aoas(maxIndex)))]);

        %% Calculate and Plot Receive Power Spectrum
        
        % Define center frequency and total bandwidth
        center_frequency = 58e9;  % Center frequency in Hz (58 GHz)
        total_bandwidth = 1e9;  % Total bandwidth in Hz (1 GHz)
        subcarrier_bandwidth = total_bandwidth / 1024;  % Bandwidth per subcarrier

        pSpectrumArray = zeros(length(tx_antennas), nFFT);  % Store power spectrum for each receive antenna

        % High-pass filter to remove DC component
        [b, a] = butter(1, 0.01, 'high'); % Butterworth filter design

        for itimes = 1:ntimes
            for rxIndex = 1:length(tx_antennas)
                td = rxtd(:, itimes, rxIndex);  % Ensure td is a vector (1024-by-1)
                td = filter(b, a, td);  % Apply high-pass filter
                fd = fft(td);  % Compute FFT (1024-by-1)
                power_spectrum = abs(fd).^2;  % Compute power spectrum (1024-by-1)
                pSpectrumArray(rxIndex, :) = pSpectrumArray(rxIndex, :) + power_spectrum';  % Sum power spectrum across receptions
            end % rxIndex
        end % itimes

        % Average the power spectrum over the number of receptions
        pSpectrumArray = pSpectrumArray / ntimes;

        % Normalize the power spectrum
        for rxIndex = 1:length(tx_antennas)
            max_power = max(pSpectrumArray(rxIndex, :));
            if max_power > 0
                pSpectrumArray(rxIndex, :) = pSpectrumArray(rxIndex, :) / max_power;
            end
        end

        % Smooth the power spectrum to reduce noise
        windowSize = 15;  % Increase window size for more smoothing
        for rxIndex = 1:length(tx_antennas)
            pSpectrumArray(rxIndex, :) = smooth(pSpectrumArray(rxIndex, :), windowSize);
        end

        % Define the frequency axis
        freq = linspace(-total_bandwidth/2, total_bandwidth/2, nFFT);

        % Convert to GHz and adjust for center frequency
        freq_GHz = center_frequency / 1e9 + freq / 1e9;

        % Plot each receive antenna's power spectrum for subcarriers -400 to 400
        figure('Position', [100, 100, figure_width, figure_height]); % Set figure size
        clf;
        hold on;
        for rxIndex = 1:length(tx_antennas)
            plot(freq_GHz(startIndex:stopIndex), mag2db(fftshift(pSpectrumArray(rxIndex, startIndex:stopIndex))), 'DisplayName', ['Ant ' num2str(tx_antennas(rxIndex))]);
        end
        hold off;

        xlabel('Frequency (GHz)', 'FontSize', FntSize);
        ylabel('Power (dB)', 'FontSize', FntSize);
        set(gca, 'FontSize', FntSize, 'YTick', common_yticks_pw);
        grid on; grid minor;
        ylim([-50 0]);
        
        % Configure legend for receive power spectrum
        lgd = legend('show');
        lgd.FontSize = 20;
        lgd.NumColumns = 2;  % Set to two columns for the receive power spectrum legend

        % Adjust legend position to fit within the plot
        lgd.Position = [0.55 - lgd.Position(3)/2, 0.25, lgd.Position(3), lgd.Position(4)];  % Center horizontally and adjust vertically

        %% Save plots
        plot_filename = ['rcv_pw_spectrum_' num2str(rx_rotation) '_' num2str(length(tx_antennas)) 'x' num2str(length(tx_antennas))];
        print(gcf, [plot_filename '.eps'], '-depsc');

    end
end
