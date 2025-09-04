%% Clear all workspace variables except for sdr objects
clc;
close all;
clear all;


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

rotation = '-30';
rotation_num = str2double(rotation); 

tx_antennas = [1, 2, 3, 4, 5, 6, 7, 8];

nFFT = 1024;
nread = nFFT;
nskip = nFFT*3;
ntimes = 100;
txfd = zeros(nFFT, 1);
scMin = -400; scMax = 400;


rx_bfm = ['rx_bfm_' num2str(rotation) '_' num2str(length(tx_antennas)) 'x' num2str(length(tx_antennas)) '.mat'];
load(rx_bfm, "rxtd");

naoa = 101;
aoas = linspace(-1.57, 1.57, naoa);
pArray = zeros(1, naoa);


for iaoa = 1:naoa
    p = 0;
    aoa = aoas(iaoa);
    for itimes = 1:ntimes
        tdbf = zeros(nFFT, 1);
        for rxIndex=1:8
            td = rxtd(:,itimes,rxIndex);
            tdbf = tdbf + td * exp(1j*rxIndex*pi*sin(aoa)); % Apply BF Vec
        end % rxIndex
        fd = fftshift(fft(tdbf));
        p = p + sum(abs(fd( nFFT/2 + 1 + scMin : nFFT/2 + 1 + scMax)));
    end %itimes
    pArray(iaoa) = p;
end % iaoa

%%  Plot
pArray = pArray / max(pArray);
figure(3); clf;
plot(rad2deg(aoas), mag2db(pArray));
xlabel('Angle of Arrival (Deg)');
ylabel('Power (dB)');
grid on; grid minor;
ylim([-12 0])
% title(sprintf('RX location %d°', rotation_num), 'FontSize', 16);

set(gca, 'FontSize', 25); % Change 12 to your desired font size

[maxValue, maxIndex] = max(pArray);
disp(['Analog Beamforming Angle: ',  num2str(rad2deg(aoas(maxIndex)))]);

plot_filename = ['beam_tracking_' num2str(rotation_num) '_' num2str(length(tx_antennas)) 'x' num2str(length(tx_antennas))];
print(gcf, [plot_filename '.eps'], '-depsc');

