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

%% Variables

character = '1';
nFFT = 1024;
tx_antennas = [1, 2, 3, 4, 5, 6, 7, 8];
name = ['gesture_CFR_' character '_' num2str(length(tx_antennas)) 'x' num2str(length(tx_antennas)) '.mat'];
load(name, "CFR_all");


CFR = squeeze(CFR_all(:, :, :, 5));

%% Plot CFR
figure()
tx_idx = 1;
for rx_idx = 1:8
    subplot(4, 2, rx_idx);
    plot(mag2db(abs(CFR(:, rx_idx, tx_idx))), 'LineWidth', 2);
    ylim([100 120]);
    xlim([0 nFFT]);
    title(strcat('CFR', num2str(tx_idx), '-', num2str(rx_idx)), 'FontSize', 20);
    set(gca, 'FontSize', 20); 
end

